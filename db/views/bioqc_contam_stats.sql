--------------------------------------------------------------------------------
-- contamination stats with meta information
--
-- List of all samples with
--   * meta information
--   * boolean flag (null | 1) indicating whether the sample is contamined. 
--------------------------------------------------------------------------------

create or replace view bioqc_contam_stats as
  with res_tissue as (
    select /*+ parallel(16) */ distinct brf.gsm
                                      , bts.tgroup
                                      , bts.tissue_set
                                      , brf.tissue
                                      , brf.organism
                                      , brf.year
                                      , brf.country
    from bioqc_res_fil brf
    join bioqc_tissue_set bts on bts.tissue = brf.tissue
  )
  select distinct br.gsm
                , br.tissue_set
                , br.tissue
                , br.organism
                , br.year
                , br.country 
                , case 
                   when cs.gsm is null
                   then null
                   else 1 
                  end as is_contam 
  from res_tissue br
  left outer join bioqc_contamined_samples cs on cs.gsm = br.gsm; 