--------------------------------------------------------------------------------
-- shortcut to detect all contamined samples. 
--
-- requires the enrichment ratio to be higher than the predefined
-- thershold for ALL expected signatuers
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
  from bioqc_res_contam brc
	group by gsm
       , tgroup
       , signature
       , signature_name
       , tissue_set
	having min(enrichment_ratio) > 6
  order by gsm, tissue_set, rk;
