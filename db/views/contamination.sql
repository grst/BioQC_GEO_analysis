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

create materialized view bioqc_res_tissue
parallel 16
build immediate
refresh force
on demand
as 
  with relevant_signatures as (
      select distinct bs.source
      from bioqc_tissue_set bts
      join bioqc_signatures bs
        on bs.id = bts.signature
  )
  select /*+ parallel(16)  */  br.* 
  from bioqc_res br
  join bioqc_signatures bs
    on bs.id = br.signature
  where bs.source in (
    --'gtex_ngs_0.85_5.gmt', 'exp.tissuemark.affy.roche.symbols.gmt'
    select  /*+ CARDINALITY(relevant_signatures, 2) */ source
    from relevant_signatures
  );
create /*+ parallel(16) */ index bioqc_res_tissue_gsm
  on bioqc_res_tissue(gsm);
create /*+ parallel(16) */ index bioqc_res_tissue_signature
  on bioqc_res_tissue(signature);



--------------------------------------------------------------------------------
-- BIOQC_RES_FIL
--
-- filtered bioqc results that can be used for analysis
--   * tissue annotated and normalized
--   * channel_count = 1 (remove obsolete studies, not considered in original
--     filtering
--   * Limit organisms to human, rat and mouse (most abundant ones) 
--
-- additionally, we add the parsed meta information year and country. 
--------------------------------------------------------------------------------

create materialized view bioqc_res_fil
parallel 16
build immediate
refresh force
on demand
as 
  select /*+ parallel(16)  */ distinct  br.gsm
                                      , br.signature
                                      , bs.name as signature_name
                                      , br.pvalue 
                                      , bnt.tissue
                                      , bg.organism_ch1 as organism
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
  from bioqc_res_tissue br
  join bioqc_gsm bg 
    on bg.gsm = br.gsm
  join bioqc_normalize_tissues bnt
    on bnt.tissue_orig = lower(bg.tissue_orig)
  join bioqc_signatures bs
    on bs.id = br.signature
  where channel_count = 1
  and bg.organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus');
  
create /*+ parallel(16) */ index bioqc_res_fil_gsm
  on bioqc_res_fil(gsm);
create /*+ parallel(16) */ index bioqc_res_fil_signature
  on bioqc_res_fil(signature);
create index bioqc_res_fil_tissue 
  on bioqc_res_fil(tissue);




--------------------------------------------------------------------------------
-- BIOQC_RES_FIL_TSET
--------------------------------------------------------------------------------

CREATE MATERIALIZED VIEW bioqc_res_fil_tset
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
  );
create /*+ parallel(16) */ index bioqc_res_fil_tset_tissue_set
  on bioqc_res_fil_tset(tissue_set);
create /*+ parallel(16) */ index bioqc_res_fil_tset_tgroup
  on bioqc_res_fil_tset(tgroup);
create /*+ parallel(16) */ index bioqc_res_fil_tset_signature
  on bioqc_res_fil_tset(signature);
  
--------------------------------------------------------------------------------
-- BIOQC_RES_CONTAM
-- 
-- for each tissue set:
--     for each sample:
--        * add the expected signatures and their pvalues
--        * calculate the enrichment ratio
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
    from bioqc_res_fil_tset brf
    
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
  -- add enrichment_ratio
  select /*+ PARALLEL(16) */ brc.*
                            , log(10, cast(
                                        exp_sig_pvalue / pvalue
                                        as binary_float))
                              as enrichment_ratio
  from bioqc_res_contam_pre brc
  order by gsm, signature; 
  
create /*+ parallel(16) */ index bioqc_res_contam_gsm
  on bioqc_res_contam(gsm);
create /*+ parallel(16) */ index bioqc_res_contam_tgroup
  on bioqc_res_contam(tgroup);
  
  
  
  
  
  
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

CREATE MATERIALIZED VIEW bioqc_contamined_samples 
parallel 16
build immediate
refresh force
on demand
as
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
  from bioqc_res_contam brc
  group by gsm
       , tgroup
       , signature
       , signature_name
       , tissue_set
  order by gsm, tissue_set, rk;
  
create /*+ parallel(16) */ index bioqc_cs_gsm
  on bioqc_contamined_samples(gsm);
create /*+ parallel(16) */ index bioqc_cs_tgroup
  on bioqc_contamined_samples(tgroup);
create /*+ parallel(16) */ index bioqc_cs_signature
  on bioqc_contamined_samples(signature);
create /*+ parallel(16) */ index bioqc_cs_min_er
  on bioqc_contamined_samples(min_enrichment_ratio);
create /*+ parallel(16) */ index bioqc_cs_rk
  on bioqc_contamined_samples(rk);
  
  
  
--------------------------------------------------------------------------------
-- BIOQC_CONTAM_STATS
-- 
-- contamination stats with meta information
--
-- List of all samples with
--   * meta information
--   * boolean flag (null | 1) indicating whether the sample is contamined. 
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
                            , cs.min_enrichment_ratio as enrichment_ratio
                            , cs.rk
                            , bts2.tgroup as destination
  from bioqc_res_fil brf
  join bioqc_tissue_set bts 
    on bts.tissue = brf.tissue
  left outer join bioqc_contamined_samples cs 
    on cs.gsm = brf.gsm 
    and cs.signature = brf.signature
    and cs.tissue_set = bts.tissue_set
  join bioqc_tissue_set bts2
    on bts2.tissue_set = bts.tissue_set
    and bts2.signature = brf.signature
  where rk = 1
--    and bts.tissue_set = 'gtex_solid'
  order by brf.gsm, cs.tissue_set, cs.rk
;

-- refresh script
BEGIN 
  DBMS_SNAPSHOT.REFRESH( '"SRSLIGHT"."BIOQC_RES_TISSUE"','C');
  DBMS_SNAPSHOT.REFRESH( '"SRSLIGHT"."BIOQC_RES_FIL_TSET"','C');
  DBMS_SNAPSHOT.REFRESH( '"SRSLIGHT"."BIOQC_RES_FIL_TSET_GTEX"','C');
  DBMS_SNAPSHOT.REFRESH( '"SRSLIGHT"."BIOQC_RES_CONTAM"','C');
  DBMS_SNAPSHOT.REFRESH( '"SRSLIGHT"."BIOQC_RES_CONTAM_GTEX"','C');
end;

