-- for migration chart with destination_tissue
with contamined_samples as (
    select gse, gsm, tissue, signature, min(enrichment_ratio) from bioqc_res_ext2 bre
    group by gse, gsm, tissue, signature
    having count(exp_sig) = (select count(*) from bioqc_tissues_signatures where tissue = bre.tissue)
) 
select cs.*, case when bts.tissue is not null then bts.tissue else 'other' end as destination_tissue from contamined_samples cs
left join bioqc_tissues_signatures bts on bts.signature = cs.signature
order by tissue 
