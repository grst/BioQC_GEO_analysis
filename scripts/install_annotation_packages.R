source("http://bioconductor.org/biocLite.R")
source("lib/db.R")

annotationPackages = dbGetQuery(mydb, "select distinct bioc_package from gpl where bioc_package is not null")
for (pkg in annotationPackages$bioc_package) {
    biocLite(sprintf("%s.db", pkg))
}
