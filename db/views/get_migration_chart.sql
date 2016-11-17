--------------------------------------------------------------------------------
-- Create view that shows migration between the tissues. 
--
--------------------------------------------------------------------------------

create view bioqc_tissue_migration as 
  with res_tissue as (
    select /*+ parallel(16) */ distinct gsm
                                      , tissue as origin
    from bioqc_res_fil 
  )
  select /*+ parallel(16) */  br.*
                            , cs.signature
                            , cs.min_enrichment_ratio
                            , case 
                                when min_enrichment_ratio is null -- sample not contamined
                                then br.origin
                                else case 
                                  when bts2.tissue is null -- sample is contamined, but with a signature that is not associated with a tissue
                                  then 'other'
                                  else bts2.tissue
                                end
                              end as destination
                            , ROW_NUMBER() over (partition by br.gsm 
                                                     order by cs.min_enrichment_ratio desc) as rk
  from res_tissue br
  left outer join bioqc_contamined_samples cs on br.gsm = cs.gsm
  left outer join bioqc_tissues_signatures bts2 on bts2.signature = cs.signature;

