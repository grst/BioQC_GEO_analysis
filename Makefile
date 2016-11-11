R=R

all: 01_create_database.html 02_select_and_get_samples.html 

clean:
	rm -fv *.html
	rm -rfv *_files

wipe: clean
	rm -rfv *_cache

01_create_database.html: 01_create_database.Rmd 
	# running create_database will mess up your db. You don't want to do that!
	pandoc -f markdown -t html -o 01_create_database.html 01_create_database.Rmd
	
02_select_and_get_samples.html: 02_select_and_get_samples.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"


