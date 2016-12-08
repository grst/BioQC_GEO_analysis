R=R
RMD_FILES= 02_select_and_get_samples.Rmd 02-2_sample_statistics.Rmd 03_make_sample_heatmaps.Rmd 04_analyse_migration.Rmd
HTML_FILES= $(patsubst %.Rmd,%.html,$(RMD_FILES))
DATA_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/gse_tissue_annot
CHUNKSUB_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/chunksub
SHELL= /bin/bash
CHUNKSUB= /pstore/home/sturmg/.local/bin/chunksub
CWD= $(shell pwd)

all: $(HTML_FILES) 01_create_database.html

.PHONY: clean
clean:
	rm -fv *.html
	rm -rfv *_files

.PHONY: wipe
wipe: clean
	rm -rfv *_cache

#################################
# Render Rmarkdown documents
#################################
01_create_database.html: 01_create_database.Rmd 
	# running create_database will mess up your db. You don't want to do that!
	# therefore, we just convert the document to html using pandoc instead of rendering. 
	pandoc -f markdown -t html -s --mathjax --number-sections -o $@ $<

$(HTML_FILES): %.html: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"




##################################
# GEO DOWNLOAD
# 
# create incremental list of files to download. 
# Then, run chunksub to download the files.
##################################
results/gse_lists/downloaded.txt: .FORCE
	find $(DATA_PATH)/geo | grep -oP "GSE(\d+)" | sort -u > $@ 

results/gse_lists/missing_download.txt: results/gse_lists/gse_tissue_annotation.txt results/gse_lists/downloaded.txt
	diff <(sort $(word 1,$^)) $(word 2,$^) | grep "^<" | grep -oP "GSE(\d+)" > $@ 

.PHONY: download_gse
download_gse: results/gse_lists/missing_download.txt
	# limit the number of concurrent jobs to 60
	$(eval CHUNKSIZE := $(shell wc -l results/gse_lists/missing_download.txt | awk '{print int($$1/60+1)}')) 
	rm -fr $(CHUNKSUB_PATH)/download_gse
	$(CHUNKSUB) -d $(CWD) -s $(CHUNKSIZE) -X y -N download_gse -j $(CHUNKSUB_PATH) "$(CWD)/scripts/geo_to_eset.R {} $(DATA_PATH)/geo" $< 




#################################
# GEO ANNOTATION
# 
# annotate expression sets with human orthologues for BioQC
#################################
results/gse_lists/annotated_esets.txt: .FORCE
	find $(DATA_PATH)/geo_annot | grep -oP "GSE(.*)\.Rdata" | sort -u > $@ 

results/gse_lists/downloaded_esets.txt: .FORCE
	find $(DATA_PATH)/geo | grep -oP "GSE(.*)\.Rdata" | sort -u > $@ 

results/gse_lists/missing_annotation.txt: results/gse_lists/downloaded_esets.txt results/gse_lists/annotated_esets.txt
	diff $^ | grep "^<" | grep -oP "GSE(.*)\.Rdata" | awk '{print "$(DATA_PATH)/geo/"$$0}' > $@

.PHONY: annotate_gse
annotate_gse: results/gse_lists/missing_annotation.txt
	$(eval CHUNKSIZE := $(shell wc -l $< | awk '{print int($$1/120+1)}'))
	rm -fr $(CHUNKSUB_PATH)/annotate_gse
	$(CHUNKSUB) -d $(CWD) -s $(CHUNKSIZE) -t /pstore/home/sturmg/.chunksub/roche_chunk.template -X y -N annotate_gse -j $(CHUNKSUB_PATH) "$(CWD)/scripts/annotate_eset.R $(DATA_PATH)/geo_annot {}" $< 



#################################
# BioQC
#
# apply BioQC to the annotated expression sets
#
# Import the bioqc_melt_all_uniq.tsv manually using Sqldeveloper. 
#################################
results/gse_lists/bioqced_esets.txt: .FORCE
	find $(DATA_PATH)/bioqc | grep -oP "GSE(.*)_bioqc_res_melt\.tab" | sed s/_bioqc_res_melt\.tab/\.Rdata/ | sort -u > $@

results/gse_lists/missing_bioqc.txt: results/gse_lists/annotated_esets.txt results/gse_lists/bioqced_esets.txt
	diff $^ | grep "^<" | grep -oP "GSE(.*)\.Rdata" | awk '{print "$(DATA_PATH)/geo_annot/"$$0}' > $@

.PHONY: run_bioqc
run_bioqc: results/gse_lists/missing_bioqc.txt results/gmt_all.gmt
	rm -fr $(CHUNKSUB_PATH)/bioqc
	$(CHUNKSUB) -d $(CWD) -s 10 -t /pstore/home/sturmg/.chunksub/roche_chunk.template -X y -N bioqc -j $(CHUNKSUB_PATH) "$(CWD)/scripts/run_bioqc.R $(DATA_PATH)/bioqc $(word 2,$^) {}" $< 

$(DATA_PATH)/bioqc_melt_all.tsv: results/gse_lists/bioqced_esets.txt 
	xargs cat < $< > $@

$(DATA_PATH)/bioqc_melt_all.uniq.tsv:  $(DATA_PATH)/bioqc_melt_all.tsv
	# unique on first two columns. For one study the results exactely identical
	# due to floating point inprecision, but identical up to 5 decimal digits (i checked)
	bash_wrapper.sh 24 "sort --parallel 24 -u -k1,2 $< > $@ 

.FORCE:
