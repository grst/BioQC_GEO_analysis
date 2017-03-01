--------------------------------------------------------------------------------
-- BIOQC_RES_TISSUE
-- 
-- Materialized view with BioQC results relevant for the BioQC GEO study. 
-- Could be replaced with a table where only relevant BioQC results
-- are imported in the first place. 
--
-- Contains all BioQC-results that have a tissue signature (ingore
-- pathway signatures).
--
-- Working on the full BioQC result table is infeasible for performance
-- reasons. 
--------------------------------------------------------------------------------

--drop materialized view bioqc_res_tissue;
--create materialized view bioqc_res_tissue
--parallel 16
--build immediate
--refresh force
--on demand
--as 
--  with relevant_signatures as (
--      select distinct bs.source
--      from bioqc_tissue_set bts
--      join bioqc_signatures bs
--        on bs.id = bts.signature
--  )
--  select /*+ parallel(16)  */  br.* 
--  from bioqc_res br
--  join bioqc_signatures bs
--    on bs.id = br.signature
--  where bs.source in (
--    --'gtex_ngs_0.85_5.gmt', 'exp.tissuemark.affy.roche.symbols.gmt'
--    select  /*+ CARDINALITY(relevant_signatures, 2) */ source
--    from relevant_signatures
--  );
--create /*+ parallel(16) */ index bioqc_res_tissue_gsm
--  on bioqc_res_tissue(gsm);
--create /*+ parallel(16) */ index bioqc_res_tissue_signature
--  on bioqc_res_tissue(signature);
  


--------------------------------------------------------------------------------
-- BIOQC_SELECTED_SAMPLES
--
-- "background" 
--------------------------------------------------------------------------------

drop materialized view bioqc_selected_samples;
create materialized view bioqc_selected_samples
parallel 16
build immediate
refresh force
on demand
as 
  select /*+ parallel(16) */ bg.gsm
                           , bg.organism_ch1 as organism
                           , bg.tissue_orig
                           , bnt.tissue
                           , cast(
                               cast(
                                 regexp_substr(submission_date, '^(\d{4})-.*', 1, 1, NULL, 1) 
                                 as varchar2(4)
                               )
                               as NUMBER(4)
                             ) as year
                           , cast(
                               TRIM(BOTH from
                                    regexp_substr(contact, 'Country:(.*?)(;.*)?$', 1, 1, NULL, 1) 
                               )
                               as varchar2(100)
                             ) as country
  from bioqc_bioqc_success bs
  join bioqc_gsm bg
    on bg.gsm = bs.gsm
  join bioqc_gse_gsm bgg
    on bgg.gsm = bs.gsm 
  join bioqc_normalize_tissues bnt
    on bnt.tissue_orig = lower(bg.tissue_orig)
  join bioqc_gse_gpl bgl
    on bgg.gse = bgl.gse
    and bg.gpl = bgl.gpl
  join bioqc_res br
    on br.gsm = bg.gsm
  where channel_count = 1
  and organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')
  -- and study_median between 3 and 9
  and ABS(study_75 - study_25) >= .5 -- IQR to ensure sufficient variance. 
  and signature = 54911 --awesome housekeepers
  and pvalue < 1e-5;
  
create /*+ parallel(16) */ index bss_gsm
  on bioqc_selected_samples(gsm); 
create /*+ parallel(16) */ index bss_tissue
  on bioqc_selected_samples(tissue);
create /*+ parallel(16) */ index bss_year
  on bioqc_selected_samples(year);
create /*+ parallel(16) */ index bss_country
  on bioqc_selected_samples(country);
  
  
