with significant_samples as (
    select distinct gsm from bioqc_res where pvalue < 1e6
) -- all samples that have at least one score >6
select * from bioqc_res r
join significant_samples ss on ss.gsm = r.gsm       -- fast filtering for siginifcant samples
join bioqc_gsm on bioqc_gsm.gsm = r.gsm             -- get tissue annotation

-- get all signatures related to tissue ("expected signatures")
join bioqc_tissues_signatures bts on bts.tissue = bioqc_gsm.tissue  

-- get pvalues of expected signatures. Match them to the respective sample.  
join bioqc_res r2 on r.gsm = r2.gsm and r.gse = r2.gse and r2.signature = bts.signature 

-- this is the threshold by which the alternate signature must exceed the expected ones. 
where r2.pvalue / r.pvalue > 1e5

-- the expected signatures are likely to meet that threshold. Filter them out. 
and r.signature not in (select signature from bioqc_tissues_signatures where tissue = bioqc_gsm.tissue) 

-- choose the tissue.
and bioqc_gsm.tissue = 'adipose'

