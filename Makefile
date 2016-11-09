R=R

all: bulk_analysis_db_prepare.html

clean:
	rm -fv *.html
	rm -rfv *_files

wipe: clean
	rm -rfv *_cache

bulk_analysis_db_prepare.html: bulk_analysis_db_prepare.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

