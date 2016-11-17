--------------------------------------------------------------------------------
-- shortcut to detect all contamined samples. 
--
-- requires the enrichment ratio to be higher than the predefined
-- thershold for ALL expected signatuers
--------------------------------------------------------------------------------
create or replace view bioqc_contamined_samples as (
	select /*+ parallel(16) */ gsm
       , tissue
       , signature
       , min(enrichment_ratio) min_enrichment_ratio 
  from bioqc_res_contam_6 bre
	group by gsm, tissue, signature
	having count(exp_sig) = (
    select count(*)
    from bioqc_tissues_signatures
    where tissue = bre.tissue
  )
)