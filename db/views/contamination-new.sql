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
-- BIOQC_SELECTED_SAMPLES
--
-- "background" 
--------------------------------------------------------------------------------

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
  join bioqc_normalize_tissues bnt
    on bnt.tissue_orig = lower(bg.tissue_orig)  
  where channel_count = 1  
  and organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus');
  
create /*+ parallel(16) */ index bss_gsm
  on bioqc_selected_samples(gsm); 
create /*+ parallel(16) */ index bss_tissue
  on bioqc_selected_samples(tissue);
create /*+ parallel(16) */ index bss_year
  on bioqc_selected_samples(year);
create /*+ parallel(16) */ index bss_country
  on bioqc_selected_samples(country);
  
  


--------------------------------------------------------------------------------
-- ## Section: Tissue Enrichment
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- BIOQC_SELECTED_SAMPLES_TSET
-- 
-- add expected signatures to selected samples
--------------------------------------------------------------------------------
drop materialized view bioqc_selected_samples_tset;
create materialized view bioqc_selected_samples_tset
parallel 16
build immediate
refresh force
on demand
as 
    select /*+ parallel(16) */   bss.gsm
                               , bts.tissue
                               , bts.tgroup
                               , bts2.signature as exp_sig
                               , bs.name as exp_sig_name
                               , case 
                                    when br.pvalue is NULL
                                    then 0.05
                                    else br.pvalue 
                                  end 
                                  as exp_sig_pvalue
                               , bts.tissue_set
    from bioqc_selected_samples bss
    -- add tgroup to sample
    join bioqc_tissue_set bts
      on bts.tissue = bss.tissue
    -- add exp_signature to tgroup
    join bioqc_tissue_set bts2
      on bts2.tgroup = bts.tgroup
      and bts2.tissue_set = bts.tissue_set
    join bioqc_signatures bs
      on bs.id = bts2.signature
    left outer join bioqc_res_tissue br 
      on br.gsm = bss.gsm 
      and br.signature = bts2.signature
    where bts.tissue_set = 'gtex_solid';  
create /*+ parallel(16) */ index bsst_gsm
  on bioqc_selected_samples_tset(gsm); 
create /*+ parallel(16) */ index bsst_tissue
  on bioqc_selected_samples_tset(tissue);
create /*+ parallel(16) */ index bsst_tgroup
  on bioqc_selected_samples_tset(tgroup);
create /*+ parallel(16) */ index bsst_signature
  on bioqc_selected_samples_tset(exp_sig);
create /*+ parallel(16) */ index bsst_tissue_set
  on bioqc_selected_samples_tset(tissue_set);


--------------------------------------------------------------------------------
-- BIOQC_RES_TSET
--
-- add tissue groups to bioqc results (pvalues)
--------------------------------------------------------------------------------

drop materialized view bioqc_res_tset;
create materialized view bioqc_res_tset
parallel 16
build immediate
refresh force
on demand
as 
  select /*+ parallel(16) */ distinct br.gsm
                           , br.signature as found_sig
                           , br.pvalue as found_sig_pvalue
                           , bs.name as found_sig_name
                           , bts.tgroup as found_tgroup
                           , bts.tissue_set
  from bioqc_res_tissue br
  join bioqc_signatures bs
    on br.signature = bs.id
  join bioqc_tissue_set bts
    on bts.signature = br.signature
  -- we can reduce the amount of data by only keeping values of selected samples
  join bioqc_selected_samples_tset bss
    on bss.tissue_set = bts.tissue_set
    and bss.gsm = br.gsm
  where bts.tissue_set = 'gtex_solid';
create /*+ parallel(16) */ index brt_gsm
  on bioqc_res_tset(gsm); 
create /*+ parallel(16) */ index brt_signature
  on bioqc_res_tset(signature);
create /*+ parallel(16) */ index brt_tgroup
  on bioqc_res_tset(tgroup);
create /*+ parallel(16) */ index brt_tissue_set
  on bioqc_res_tset(tissue_set);
  
--------------------------------------------------------------------------------


drop materialized view bioqc_tissue_enrichment;
create materialized view bioqc_tissue_enrichment
parallel 16
build immediate
refresh force
on demand
as 
  with bioqc_tissue_enrichment_1 as (
    select /*+ parallel(16) */  bss.gsm
                               , bss.tissue_set
                               , bss.tissue
                               , bss.tgroup
                               , bss.exp_sig
                               , bss.exp_sig_name
                               , bss.exp_sig_pvalue
                               , br.found_sig
                               , br.found_sig_pvalue
                               , br.found_sig_name
                               , br.found_tgroup
    from bioqc_selected_samples_tset bss
    left outer join bioqc_res_tset br
      on br.gsm = bss.gsm 
      and br.tissue_set = bss.tissue_set
    -- exclude expected signatures from found enriched signatures
    -- (e.g. tissue is colon, jejunum is found enriched but we don't see 
    -- that as contamination
    where br.found_sig not in (
      select signature
      from bioqc_tissue_set bts
      where bts.tissue_set = bss.tissue_set
        and bts.tgroup = bss.tgroup
    )
  ) 
  select /*+ PARALLEL(16) */ bre.*
                            , log(10, cast(
                                        exp_sig_pvalue / found_sig_pvalue
                                        as binary_float))
                              as enrichment_ratio
  from bioqc_tissue_enrichment_1 bre;
 
create /*+ parallel(16) */ index bte_gsm
  on bioqc_tissue_enrichment(gsm); 
create /*+ parallel(16) */ index bte_tissue_set
  on bioqc_tissue_enrichment(tissue_set);
create /*+ parallel(16) */ index bte_tgroup
  on bioqc_tissue_enrichment(tgroup);
create /*+ parallel(16) */ index bte_exp_sig
  on bioqc_tissue_enrichment(exp_sig);
create /*+ parallel(16) */ index bte_found_sig
  on bioqc_tissue_enrichment(found_sig);
create /*+ parallel(16) */ index bte_found_sig_name
  on bioqc_tissue_enrichment(found_sig_name);
create /*+ parallel(16) */ index bte_found_tgroup
  on bioqc_tissue_enrichment(found_tgroup);
create /*+ parallel(16) */ index bte_enrichment_ratio
  on bioqc_tissue_enrichment(enrichment_ratio);
  
  
  
