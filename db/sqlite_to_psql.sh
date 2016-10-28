#!/bin/bash 

#########
# Bioqc Database Creation
#
# In this file I describe step by step
# how I created the BioQC GEO databse. 
echo "This bash script is not meant to be executed as a script,
it rather describes line by line how the bioqc database was created." 
read -n1 -r -p "Press any key to continue.." key
########


# first retrieve the database schema: 
sqlite3 GEOmetadb.sqlite .schema > geometadb_schema.sql 

# create schema in psql
# 
# It might be the case, that we have to fix that sql file manually first. 
psql -f geometadb_schema.sql 

# export the tables 
mkdir -p tables
for table in $(cat tables.txt); do 
    sqlite3 GEOmetadb.sqlite -header -csv -separator ',' "select * from ${table};" > tables/${table}.csv;
done

# remove invalid utf8 characters
for file in tables/*.csv; do 
    iconv -c -f utf-8 -t utf-8 $file > $file.utf8.csv; 
done 

# import the tables into postgres
for table in $(cat tables.txt); do
    echo ${table};
    psql -c "\copy ${table} from tables/${table}.csv.utf8.csv with csv header encoding 'UTF8' delimiter as ','"; 
done

# for some reason it did not work to import the gsm table as csv file
# (random lines were missing, apparently due to a kernel bug)
# 
# I therefore imported the gsm table by dumping the sqlite sql code
# and executing the insert statements manually. 
# Actually I had to run it twice to get all entries. 
sqlite3 GEOmetadb.sqlite ".dump gsm" > tables/gsm.sql 
iconv -c -f utf-8 -t utf-8 tables/gsm.sql > tables/gsm.sql.utf8.sql
psql -f <(grep "^INSERT" tables/gsm.sql.utf8.sql)

# apply the schema modifications.
# 
# add foreign keys et cetera. 
psql -f update_geometabase.sql



