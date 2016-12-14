--------------------------------------------------------------------------------
-- Create view that shows migration between the tissues. 
--
--------------------------------------------------------------------------------

create or replace view bioqc_tissue_migration as 
  with res_tissue as (
    select /*+ parallel(16) */ distinct gsm
    from bioqc_res_fil 
  )
  select /*+ parallel(16) */  br.gsm
                            , br.signature
                            , br.signature_name
                            , br.pvalue
                            , br.organism
                            , br.year
                            , br.country
                            , bts.tgroup as origin 
                            , bts.tissue_set
                            , cs.min_enrichment_ratio
                            , cs.rk
                            , case 
                                when min_enrichment_ratio is null -- sample not contamined
                                then bts.tgroup
                                else case 
                                  when bts2.tgroup is null -- sample is contamined, but with a signature that is not associated with a tissue
                                  then 'other'
                                  else bts2.tgroup
                                end
                              end as destination
  from bioqc_res_fil br
  join bioqc_tissue_set bts 
    on bts.tissue = br.tissue
  left outer join bioqc_contamined_samples cs 
    on cs.gsm = br.gsm 
    and cs.signature = br.signature
    and cs.tissue_set = bts.tissue_set
  join bioqc_tissue_set bts2
    on bts2.tissue_set = bts.tissue_set
    and bts2.signature = br.signature
  where (rk is null or rk = 1)
     -- and bts.tissue_set = 'gtex_solid'
  order by br.gsm, cs.tissue_set, cs.rk
;

