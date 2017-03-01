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
                               , bts.tissue_set
                               , bts.tgroup
                               , bts2.signature as exp_sig
                               , bs.name as exp_sig_name
                               , case 
                                    when br.pvalue is NULL
                                    then 0.1
                                    else br.pvalue 
                                  end 
                                  as exp_sig_pvalue
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
    left outer join bioqc_res br 
      on br.gsm = bss.gsm 
      and br.signature = bts2.signature;
--    where bts.tissue_set = 'gtex_solid';  
create /*+ parallel(16) */ index bioqc_sst_gsm
  on bioqc_selected_samples_tset(gsm); 
create /*+ parallel(16) */ index bioqc_sst_tissue
  on bioqc_selected_samples_tset(tissue);
create /*+ parallel(16) */ index bioqc_sst_tgroup
  on bioqc_selected_samples_tset(tgroup);
create /*+ parallel(16) */ index bioqc_sst_signature
  on bioqc_selected_samples_tset(exp_sig);
create /*+ parallel(16) */ index bioqc_sst_tissue_set
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
  from bioqc_res br
  -- we can reduce the amount of data by only keeping values of selected samples
  join bioqc_selected_samples bss
    on bss.gsm = br.gsm
  join bioqc_signatures bs
    on br.signature = bs.id
  join bioqc_tissue_set bts
    on bts.signature = br.signature;
--  where bts.tissue_set = 'gtex_solid';
create /*+ parallel(16) */ index bioqc_rt_gsm
  on bioqc_res_tset(gsm); 
create /*+ parallel(16) */ index bioqc_rt_found_sig
  on bioqc_res_tset(found_sig);
create /*+ parallel(16) */ index bioqc_rt_found_tgroup
  on bioqc_res_tset(found_tgroup);
create /*+ parallel(16) */ index bioqc_rt_tissue_set
  on bioqc_res_tset(tissue_set);
  


--------------------------------------------------------------------------------
-- BIOQC_TISSUE_ENRICHMENT
-- 
-- put samples and result together. 
-- calculate enrichment score for each signature
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
    where br.found_sig is null 
      or br.found_sig not in (
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
 
create /*+ parallel(16) */ index bioqc_te_gsm
  on bioqc_tissue_enrichment(gsm); 
create /*+ parallel(16) */ index bioqc_te_tissue_set
  on bioqc_tissue_enrichment(tissue_set);
create /*+ parallel(16) */ index bioqc_te_tgroup
  on bioqc_tissue_enrichment(tgroup);
create /*+ parallel(16) */ index bioqc_te_exp_sig
  on bioqc_tissue_enrichment(exp_sig);
create /*+ parallel(16) */ index bioqc_te_found_sig
  on bioqc_tissue_enrichment(found_sig);
create /*+ parallel(16) */ index bioqc_te_found_sig_name
  on bioqc_tissue_enrichment(found_sig_name);
create /*+ parallel(16) */ index bioqc_te_found_tgroup
  on bioqc_tissue_enrichment(found_tgroup);
create /*+ parallel(16) */ index bioqc_te_enrichment_ratio
  on bioqc_tissue_enrichment(enrichment_ratio);



--------------------------------------------------------------------------------
-- BIOQC_TISSUE_ENRICHMENT2
--
-- combine expected signature by taking the minimal enrichment ratio
-- for each expected signature. 
--
-- Example:
-- GSM    tissue    tgroup    expected  found         enrichment_ratio
-- GSM888 jejunum   intestine colon     liver_fetal   12
-- GSM888 jejunum   intestine colon     liver         8
-- GSM888 jejunum   intestine jejunum   liver_fetal   5
-- GSM888 jejunum   intestine jejunum   liver         4
--
-- will be combined into
-- GSM888           intestine           liver_fetal   5
-- GSM888           intestine           liver         4
--
-- if multiple infiltrating tissues are found, a rank is calculated. 
--------------------------------------------------------------------------------

drop materialized view bioqc_tissue_enrichment2;
CREATE MATERIALIZED VIEW bioqc_tissue_enrichment2 
parallel 16
build immediate
refresh force
on demand
as  
 select /*+ parallel(16) */ bre.gsm
                          , bre.tissue_set
                          , bre.tgroup
                          , bre.found_sig
                          , bre.found_sig_name
                          , min(enrichment_ratio) as min_enrichment_ratio
                          , ROW_NUMBER() over (
                             partition by gsm, tgroup, tissue_set 
                             order by min(enrichment_ratio) desc)
                             as rk
    
  from bioqc_tissue_enrichment bre
  group by bre.gsm
         , bre.tissue_set
         , bre.tgroup
         , bre.found_sig
         , bre.found_sig_name
  order by gsm, tissue_set, rk;
  
create /*+ parallel(16) */ index bioqc_te2_gsm
  on bioqc_tissue_enrichment2(gsm);
create /*+ parallel(16) */ index bioqc_te2_tgroup
  on bioqc_tissue_enrichment2(tgroup);
create /*+ parallel(16) */ index bioqc_te2_found_sig
  on bioqc_tissue_enrichment2(found_sig);
create /*+ parallel(16) */ index bioqc_te2_min_er
  on bioqc_tissue_enrichment2(min_enrichment_ratio);
create /*+ parallel(16) */ index bioqc_te2_rk
  on bioqc_tissue_enrichment2(rk);
  
-------------------------------------------------------------------------------
-- BIOQC_CONTAM_STATS
-- 
-- contamination stats with meta information
-- 
--------------------------------------------------------------------------------
create or replace view bioqc_contam_stats
as
  select /*+ parallel(16) */ distinct bss.gsm
                , bte.tissue_set
                , bss.tissue
                , bte.tgroup
                , bte.found_sig
                , bte.found_sig_name
                , bte.min_enrichment_ratio as enrichment_ratio
                , bss.organism
                , bss.year
                , bss.country
  from bioqc_selected_samples bss
  join bioqc_tissue_enrichment2 bte
    on bte.gsm = bss.gsm
  where rk = 1;
  
--------------------------------------------------------------------------------
-- BIOQC_TISSUE_MIGRATION
--
-- Create view that shows migration between the tissues. 
--------------------------------------------------------------------------------
create or replace view bioqc_tissue_migration
as
  select /*+ parallel(16) */  distinct bte.gsm
                          , bte.tissue_set
                          , bte.tgroup as origin
                          , bte.found_sig
                          , bte.found_sig_name
                          , bte.min_enrichment_ratio as enrichment_ratio
                          , bte.rk 
                          , case when bts.tgroup is null
                              then bte.tgroup
                              else bts.tgroup
                            end as destination
  from bioqc_tissue_enrichment2 bte
  left outer join bioqc_tissue_set bts
    on bts.signature = bte.found_sig
    and bts.tissue_set = bte.tissue_set;
  
