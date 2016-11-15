R=R
RMD_FILES= 02_select_and_get_samples.Rmd 03_analyse_samples.Rmd
HTML_FILES= $(patsubst %.Rmd,%.html,$(RMD_FILES))

all: $(HTML_FILES) 01_create_database.html

clean:
	rm -fv *.html
	rm -rfv *_files

wipe: clean
	rm -rfv *_cache

01_create_database.html: 01_create_database.Rmd 
	# running create_database will mess up your db. You don't want to do that!
	pandoc -f markdown -t html -s --mathjax --number-sections -o $@ $<
	
$(HTML_FILES): %.html: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

