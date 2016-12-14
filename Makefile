R=R
RMD_FILES= 01_create_database.Rmd 02_select_and_get_samples.Rmd 02-2_sample_processing.Rmd 03_make_sample_heatmaps.Rmd 04_analyse_migration.Rmd
PREVIEW_FILES = $(patsubst %,%.preview,$(RMD_FILES))
DATA_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/gse_tissue_annot
CHUNKSUB_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/chunksub
SHELL= /bin/bash
CHUNKSUB= /pstore/home/sturmg/.local/bin/chunksub
CWD= $(shell pwd)


#################################
# Render Rmarkdown documents
#################################

.PHONY: book
book: $(RMD_FILES)
	Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"

.PHONY: upload-book
upload-book: book
	cd gh-pages && cp -R ../_book/* ./ && git add --all * && git commit --allow-empty -m "update docs" && git push origin gh-pages

# render a chapter only by calling `make chapter1.Rmd.preview`
.PHONY: $(PREVIEW_FILES)
$(PREVIEW_FILES): %.Rmd.preview: %.Rmd
	Rscript -e "bookdown::preview_chapter('$<', 'bookdown::gitbook')"

.PHONY: clean
clean:
	rm -fv *.html
	rm -fv *.CommonMark
	rm -rfv *_files
	rm -rfv *_book 
	rm -fv _main*

.PHONY: wipe
wipe: clean
	rm -rfv *_cache



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
	find $(DATA_PATH)/bioqc | grep -oP "GSE(.*)_bioqc_res_melt\.tab" | sort -u > $@

results/gse_lists/missing_bioqc.txt: results/gse_lists/annotated_esets.txt results/gse_lists/bioqced_esets.txt
	diff $(word 1,$^) <(sed s/_bioqc_res_melt\.tab/\.Rdata/ $(word 2,$^)) | grep "^<" | grep -oP "GSE(.*)\.Rdata" | awk '{print "$(DATA_PATH)/geo_annot/"$$0}' > $@

.PHONY: run_bioqc
run_bioqc: results/gse_lists/missing_bioqc.txt results/gmt_all.gmt
	rm -fr $(CHUNKSUB_PATH)/bioqc
	$(CHUNKSUB) -d $(CWD) -s 10 -t /pstore/home/sturmg/.chunksub/roche_chunk.template -X y -N bioqc -j $(CHUNKSUB_PATH) "$(CWD)/scripts/run_bioqc.R $(DATA_PATH)/bioqc $(word 2,$^) {}" $< 

$(DATA_PATH)/bioqc_melt_all.tsv: results/gse_lists/bioqced_esets.txt 
	awk '{print "$(DATA_PATH)/bioqc/"$$0}' < $< | xargs cat > $@

$(DATA_PATH)/bioqc_melt_all.uniq.tsv:  $(DATA_PATH)/bioqc_melt_all.tsv
	# unique on first two columns. For one study the results exactely identical
	# due to floating point inprecision, but identical up to 5 decimal digits (i checked)
	bash_wrapper.sh 24 "sort --parallel 24 -u -k1,2 $< > $@ "

.FORCE:
