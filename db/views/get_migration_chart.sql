--------------------------------------------------------------------------------
-- Create view that shows migration between the tissues. 
--
--------------------------------------------------------------------------------

create or replace view bioqc_tissue_migration as 
  with res_tissue as (
    select /*+ parallel(16) */ distinct gsm
    from bioqc_res_fil 
  )
  select /*+ parallel(16) */  br.*
                            , cs.tissue_set 
                            , cs.tgroup as origin
                            , cs.signature
                            , cs.signature_name
                            , cs.min_enrichment_ratio
                            , cs.rk
                            , case 
                                when min_enrichment_ratio is null -- sample not contamined
                                then cs.tgroup
                                else case 
                                  when bts2.tgroup is null -- sample is contamined, but with a signature that is not associated with a tissue
                                  then 'other'
                                  else bts2.tgroup
                                end
                              end as destination
  from res_tissue br
  left outer join bioqc_contamined_samples cs 
    on br.gsm = cs.gsm
  left outer join bioqc_tissue_set bts2 
    on bts2.signature = cs.signature
    and bts2.tissue_set = cs.tissue_set
--  where cs.tissue_set = 'gtex_solid'
  order by br.gsm, cs.tissue_set, cs.rk
;

