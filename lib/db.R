require("RJDBC")
source("lib/db/.db_creds.R")
.jinit()
options(java.parameters = "-Xmx4g" )
drv <- RJDBC::JDBC("org.postgresql.Driver", "lib/db/postgresql-9.4.1211.jre7.jar")
mydb <- dbConnect(drv, sprintf("jdbc:postgresql://%s:5432/%s", pg_server, pg_dbname), user = pg_user, password = pg_pass )

#' Wapper for \code{\link{dbWriteTable}}.
#'
#' Calls dbWriteTable with append=TRUE and overwrite=FALSE as defaults. 
#' 
#' @param df
#' @param table
dbAppendDf = function(table, df) {
  dbWriteTable(mydb, name=table, value=df, append=TRUE, overwrite=FALSE)
}