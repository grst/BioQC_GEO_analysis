stopifnot(suppressPackageStartupMessages(require(ribiosUtils)))
source("lib/db/.db_creds.R")
mydb = ribiosCon(db=pg_dbname, user=pg_user, password=pg_pass, forceJDBC=TRUE)

#' Wapper for \code{\link{dbWriteTable}}.
#'
#' Calls dbWriteTable with append=TRUE and overwrite=FALSE as defaults. 
#' 
#' @param df
#' @param table
dbAppendDf = function(table, df) {
  #tmp = tempfile()
  #write.csv(df, file=tmp, fileEncoding='utf-8')
 # system(sprintf(
 #  'PGPASSFILE=/homebasel/biocomp/sturmg/.pgpass /apps64/postgresql-9.2.2/bin/psql -U sturmg -c "\\copy %s from %s with csv header encoding \'UTF8\' delimiter as \',\'"', table, tmp))
  dbWriteTable(mydb, name=table, value=df, append=TRUE, overwrite=FALSE)
}