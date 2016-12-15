CREATE MATERIALIZED VIEW bioqc_res_fil_tset_gtex
parallel 16
build immediate
refresh force
on demand
as
 -- distinct, because the 'expected tissue' is omitted 
 -- which would result in duplicated rows
 select /*+ parallel(16) */ distinct  brf.gsm
                                    , brf.tissue
                                    , bts.tgroup
                                    , bts.tissue_set 
                                    , brf.signature
                                    , brf.signature_name
                                    , brf.pvalue
  from bioqc_res_fil brf
  join bioqc_tissue_set bts 
    on bts.tissue = brf.tissue
  -- make sure only the signature from the same signature_set are included. 
  where brf.signature in (
    select bts2.signature
    from bioqc_tissue_set bts2
    where bts2.tissue_set = bts.tissue_set
  )
  and tissue_set = 'gtex_solid';
create /*+ parallel(16) */ index brft_gtex_tissue_set
  on bioqc_res_fil_tset_gtex(tissue_set);
create /*+ parallel(16) */ index brft_gx_tgroup
  on bioqc_res_fil_tset_gtex(tgroup);
create /*+ parallel(16) */ index brft_gx_signature
  on bioqc_res_fil_tset_gtex(signature);
  
  
CREATE MATERIALIZED VIEW bioqc_res_contam_gtex 
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
    join bioqc_signatures bs
      on bs.id = bts.signature
  ),
  bioqc_res_contam_pre as (
    -- distinct because the tgroup in expected_signatures has multiple
    -- tissues which are omitted and therefore lead to duplicate rows. 
    select /*+ PARALLEL(16) */ distinct brf.gsm
                                      , brf.tissue
                                      , brf.tgroup
                                      , brf.tissue_set
                                      , brf.signature
                                      , brf.signature_name
                                      , brf.pvalue                            
                                      , es.signature as exp_sig
                                      , es.signature_name as exp_sig_name
                                      -- if no pvalue existant for expected signature, assume that it is available at the cutoff value 0.05
                                      , case 
                                          when brf2.pvalue is NULL
                                          then 0.05
                                          else brf2.pvalue 
                                        end 
                                        as exp_sig_pvalue
    from bioqc_res_fil_tset_gtex brf
    
    -- get all signatures related to tissue groups ("expected signatures")
    join expected_signatures es 
      on es.tgroup = brf.tgroup
      and es.tissue_set = brf.tissue_set
    
    -- get pvalues of expected signatures. Match them to the respective sample.  
    left outer join bioqc_res_fil brf2 
      on brf2.gsm = brf.gsm 
      and brf2.signature = es.signature 
      
    -- exclude expected signatures from found enriched signatures
    -- (e.g. tissue is colon, jejunum is found enriched but we don't see 
    -- that as contamination
    where brf.signature not in (
      select signature
      from expected_signatures es2
      where es2.tgroup = brf.tgroup
      and es2.tissue_set = brf.tissue_set
    )
  )
  select /*+ PARALLEL(16) */ brc.*
                            , log(10, cast(
                                        exp_sig_pvalue / pvalue
                                        as binary_float))
                              as enrichment_ratio
  from bioqc_res_contam_pre brc
  order by gsm, signature; 
  
create /*+ parallel(16) */ index bioqc_rescontamg_gsm
  on bioqc_res_contam_gtex(gsm);
create /*+ parallel(16) */ index bioqc_rescontamg_tgroup
  on bioqc_res_contam_gtex(tgroup);
  



--------------------------------------------------------------------------------
-- BIOQC_CONTAMINED_SAMPLES
--
-- View containing all samples we defined as contamined 
-- (enrichment_ratio > 6)
--
-- requires the enrichment ratio to be higher than the predefined
-- thershold for ALL expected signatuers
--
-- The only point of group by is creating the rank
--------------------------------------------------------------------------------
create or replace view bioqc_contamined_samples as 
	select /*+ parallel(16) */ gsm
       , tgroup
       , signature
       , signature_name
       , tissue_set
       , min(enrichment_ratio) min_enrichment_ratio 
       , ROW_NUMBER() over (
           partition by gsm, tgroup, tissue_set 
           order by min(enrichment_ratio) desc)
           as rk
  from bioqc_res_contam_gtex brc
	group by gsm
       , tgroup
       , signature
       , signature_name
       , tissue_set
  order by gsm, tissue_set, rk;






--------------------------------------------------------------------------------
-- BIOQC_CONTAM_STATS
-- 
-- contamination stats with meta information
--
-- List of all samples with
--   * meta information
--   * enrichment ratio, if available
--
-- if multiple signatures show up, the maximal enrichment ratio is taken. 
--------------------------------------------------------------------------------

create or replace view bioqc_contam_stats as
  select distinct brf.gsm
                , bts.tissue_set
                , bts.tissue
                , brf.organism
                , brf.year
                , brf.country 
                , cs.min_enrichment_ratio as enrichment_ratio
                , cs.signature
                , cs.signature_name
  from bioqc_res_fil brf
  join bioqc_tissue_set bts 
    on bts.tissue = brf.tissue
  left outer join bioqc_contamined_samples cs
    on cs.gsm = brf.gsm
    and cs.tissue_set = bts.tissue_set
  where rk = 1;
  
  
--------------------------------------------------------------------------------
-- BIOQC_TISSUE_MIGRATION
--
-- Create view that shows migration between the tissues. 
--------------------------------------------------------------------------------

create or replace view bioqc_tissue_migration as 
  with res_tissue as (
    select /*+ parallel(16) */ distinct gsm
    from bioqc_res_fil 
  )
  select /*+ parallel(16) */  brf.gsm
                            , brf.signature
                            , brf.signature_name
                            , brf.pvalue
                            , brf.organism
                            , brf.year
                            , brf.country
                            , bts.tgroup as origin 
                            , bts.tissue_set
                            , cs.min_enrichment_ratio
                            , cs.rk
                            , bts2.tgroup as destination
  from bioqc_res_fil brf
  join bioqc_tissue_set bts 
    on bts.tissue = br.tissue
  left outer join bioqc_contamined_samples cs 
    on cs.gsm = br.gsm 
    and cs.signature = br.signature
    and cs.tissue_set = bts.tissue_set
  join bioqc_tissue_set bts2
    on bts2.tissue_set = bts.tissue_set
    and bts2.signature = br.signature
  where rk = 1
  order by br.gsm, cs.tissue_set, cs.rk
;

 