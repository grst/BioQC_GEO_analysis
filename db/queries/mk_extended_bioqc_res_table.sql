-------------------------------------------------------------------------------------------
-- create table (to be seen as materialized view, which does not exist in postgres 9.2)
-- which contains the enrichment ratio for all siginificant samples and 
-- all expected signatures for each tissue. 
--
-- this is *very* slow as it creates a *huge* table with >80M row. 
-------------------------------------------------------------------------------------------

CREATE TABLE bioqc_res_ext as
with significant_samples as (
    select distinct gse, gsm from bioqc_res where pvalue < 1e-6
) -- all samples that have at least one score >6

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

where (r.gse, r.gsm) in (select * from significant_samples)    	-- filtering for significant samples

-- the expected signatures are likely to meet that threshold. Filter them out. 
and r.signature not in (select signature from bioqc_tissues_signatures where tissue = bioqc_gsm.tissue) 



