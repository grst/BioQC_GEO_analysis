--------------------------------------------------------------------------------
-- contamination stats with meta information
--
-- List of all samples with
--   * meta information
--   * boolean flag (null | 1) indicating whether the sample is contamined. 
--------------------------------------------------------------------------------

create view bioqc_contam_stats as
  select distinct brf.gsm
                , brf.tissue
                , brf.organism
                , brf.year
                , brf.country 
                , case 
                   when cs.gsm is null
                   then null
                   else 1 
                  end as is_contam 
  from bioqc_res_fil brf
  left outer join bioqc_contamined_samples cs on cs.gsm = brf.gsm 