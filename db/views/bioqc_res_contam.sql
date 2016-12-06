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

CREATE MATERIALIZED VIEW bioqc_res_contam 
parallel 16
build immediate
refresh force
on demand
as
  with 
  expected_signatures as (
    select  bts.* 
          , bs.name as signature_name
    from bioqc_tissue_set bts
    join bioqc_signatures bs on bs.id = bts.signature
  ), 
  bioqc_res_fil_tset as (
    select  brf.gsm
          , brf.tissue
          , bts.tgroup
          , bts.tissue_set 
          , brf.signature
          , brf.signature_name
          , brf.pvalue
    from bioqc_res_fil brf
    join bioqc_tissue_set bts 
      on bts.tissue = brf.tissue
    -- make sure only the signature from the same signature set are included. 
    where brf.signature in (
      select bts2.signature
      from bioqc_tissue_set bts2
      where bts2.tissue_set = bts.tissue_set
    )
  )
  select /*+ PARALLEL(16) */  distinct -- distinct, because tissue -> tissue_group -> tissue is not unique. 
                              brf.gsm
                            , brf.tissue
                            , brf.tgroup
                            , brf.tissue_set
                            , brf.signature
                            , brf.signature_name
                            , brf.pvalue                            
                            , es.signature as exp_sig
                            , es.signature_name as exp_sig_name
                            , brf2.pvalue as exp_sig_pvalue
                            , log(10, cast(brf2.pvalue / brf.pvalue as binary_float)) as enrichment_ratio 
  from bioqc_res_fil_tset brf
  
  -- get all signatures related to tissue groups ("expected signatures")
  join expected_signatures es 
    on es.tgroup = brf.tgroup
    and es.tissue_set = brf.tissue_set
  
  -- get pvalues of expected signatures. Match them to the respective sample.  
  join bioqc_res_fil_tset brf2 
    on brf2.gsm = brf.gsm 
    and brf2.signature = es.signature 
    
  -- exclude expected signatures (e.g. jejunum having a high enrichment ratio
  -- in relation to colon)
  where brf.signature not in (
    select signature
    from expected_signatures es2
    where es2.tgroup = brf.tgroup
    and es2.tissue_set = brf.tissue_set
  )
  
--  where 1 = 1
  order by gsm, brf.signature;