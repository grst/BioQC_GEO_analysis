R=R
RMD_FILES= 02_select_and_get_samples.Rmd 02-2_sample_statistics.Rmd 03_make_sample_heatmaps.Rmd 04_analyse_migration.Rmd
HTML_FILES= $(patsubst %.Rmd,%.html,$(RMD_FILES))
DATA_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/gse_tissue_annot
CHUNKSUB_PATH= /pstore/data/biocomp/users/sturmg/BioQC_GEO_analysis/chunksub
SHELL= /bin/bash
CHUNKSUB= /pstore/home/sturmg/.local/bin/chunksub
CWD= $(shell pwd)

all: $(HTML_FILES) 01_create_database.html

clean:
	rm -fv *.html
	rm -rfv *_files

wipe: clean
	rm -rfv *_cache

#################################
# Render Rmarkdown documents
# ###############################
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
results/gse_lists/downloaded.txt: 
	find $(DATA_PATH)/geo | grep -oP "GSE(\d+)" | sort -u > results/gse_lists/downloaded.txt

results/gse_lists/missing_download.txt: results/gse_lists/gse_tissue_annotation.txt results/gse_lists/downloaded.txt
	diff <(sort results/gse_lists/gse_tissue_annotation.txt) results/gse_lists/downloaded.txt | grep "^<" | grep -oP "GSE(\d+)" > results/gse_lists/missing_download.txt 

download_gse: results/gse_lists/missing_download.txt
	# limit the number of concurrent jobs to 60
	$(eval CHUNKSIZE := $(shell wc -l results/gse_lists/missing_download.txt | awk '{print int($$1/60+1)}')) 
	$(CHUNKSUB) -d $(CWD) -s $(CHUNKSIZE) -X y -N download_gse -j $(CHUNKSUB_PATH) "$(CWD)/scripts/geo_to_eset.R {} $(DATA_PATH)/geo" results/gse_lists/missing_download.txt

# Annotate ExpressionSets with orthologous gene symbols for BioQC

# run BioQC
bioqc_melet_list.txt: geo 
	find `pwd`/geo -type f > eset_list.txt 

eset_annot_list.txt: geo 
	find `pwd`/geo_annot -type f > eset_annot_list.txt 
	.tsv: bioqc_all
	cat bioqc_all/*_melt.tab > bioqc_melt_all.tsv
	# unique on first two columns. For one study the results exactely identical
	# due to floating point inprecision, but identical up to 5 decimal digits (i checked)
	./bash_wrapper.sh 24 "sort --parallel 24 -u -k1,2 bioqc_melt_all.tsv > bioqc_melt_all.uniq.tsv"


