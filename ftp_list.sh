#!/bin/bash
ml load gnutls
/pstore/home/sturmg/bin/lftp -c "open ftp.ncbi.nlm.nih.gov/geo/platforms/${1} && find && exit" > $1.out.txt
