--------------------------------------------------------------------------------
-- create table (to be seen as materialized view, which does not exist in postgres 9.2)
-- which contains all siginificant samples and enrichment_ratios for all
-- expected signatures, which have an enrichment ratio higher than a certain threshold. 
--
-- significantly faster than creating the whole table without enrichment-ratio cutoff
-- with enrichment-ratio > 6 -> 300k instead of 80M rows. 
--------------------------------------------------------------------------------

CREATE TABLE bioqc_res_ext2 as
select 	r.*, 
		bioqc_gsm.tissue as tissue, 
		bts.signature as exp_sig, 
		r2.pvalue as exp_sig_pval,
		log(10, cast(r2.pvalue / r.pvalue as numeric)) as enrichment_ratio 
from bioqc_res r
join bioqc_gsm on bioqc_gsm.gsm = r.gsm             					-- get tissue annotation

-- get all signatures related to tissue ("expected signatures")
join bioqc_tissues_signatures bts on bts.tissue = bioqc_gsm.tissue  

-- get pvalues of expected signatures. Match them to the respective sample.  
join bioqc_res r2 on r.gsm = r2.gsm and r.gse = r2.gse and r2.signature = bts.signature 

-- the expected signatures are likely to meet that threshold. Filter them out. 
where r.signature not in (select signature from bioqc_tissues_signatures where tissue = bioqc_gsm.tissue) 
and r2.pvalue / r.pvalue > 1e6


