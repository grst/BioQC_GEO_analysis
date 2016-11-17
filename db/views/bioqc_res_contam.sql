--------------------------------------------------------------------------------
-- create materialized view 
-- which contains all siginificant samples and enrichment_ratios for all
-- expected signatures, which have an enrichment ratio higher 
-- than a certain threshold. 
--
-- significantly faster than creating the whole table without 
-- enrichment-ratio cutoff with enrichment-ratio > 6 
-- (-> 300k instead of 80M rows)
--------------------------------------------------------------------------------

CREATE MATERIALIZED VIEW bioqc_res_contam_6 -- log fold > 6
parallel 16
build immediate
refresh force
on demand
as
  select /*+ PARALLEL(16) */  brf.gsm
                            , brf.signature
                            , brf.pvalue
                            , brf.tissue
                            , bts.signature as exp_sig
                            , brf2.pvalue as exp_sig_pval
                            , log(10, cast(brf2.pvalue / brf.pvalue as binary_float)) as enrichment_ratio 
  from bioqc_res_fil brf
  
  -- get all signatures related to tissue ("expected signatures")
  join bioqc_tissues_signatures bts on bts.tissue = brf.tissue 
  
  -- get pvalues of expected signatures. Match them to the respective sample.  
  join bioqc_res brf2 on brf.gsm = brf2.gsm and brf2.signature = bts.signature 
  
  -- the expected signatures are likely to meet that threshold. Filter them out. 
  where brf.signature not in (
    select signature
    from bioqc_tissues_signatures
    where tissue = brf.tissue
  ) 
  
  --- TODO
  --- select for significant samples here, if enrichment-ratio < 6
  ---
  --- filter for enrichment ratio > log6 fold
  and brf2.pvalue / brf.pvalue > 1e6;
  
create /*+ parallel(32) */ index idx_bioqc_res_contam_6_gsm on bioqc_res_contam_6(gsm);
create /*+ parallel(32) */ index bioqc_res_contam_6_signature on bioqc_res_contam_6(signature);
