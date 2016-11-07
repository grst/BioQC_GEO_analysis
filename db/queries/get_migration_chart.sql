-- for migration chart with destination_tissue
with contamined_samples as (
    select gse, gsm, tissue, signature, min(enrichment_ratio), count(exp_sig) as count_exp from bioqc_res_ext2 bre
    group by gse, gsm, tissue, signature
),
migration as (
	select cs.*,
	case when count_exp = (
		select count(*) from bioqc_tissues_signatures bts2
		where bts2.tissue = cs.tissue)
	then case when bts.tissue is not null
		then bts.tissue
		else 'other'
		end
	else cs.tissue
	end as destination_tissue from contamined_samples cs
	left join bioqc_tissues_signatures bts on bts.signature = cs.signature
)
select tissue, destination_tissue, count(gsm) as "count"
from migration
group by tissue, destination_tissue


-- for migration chart with destination_tissue
-- with contamined_samples as (
--    select gse, gsm, tissue, signature, min(enrichment_ratio) from bioqc_res_ext2 bre
--    group by gse, gsm, tissue, signature
--    having count(exp_sig) = (select count(*) from bioqc_tissues_signatures where tissue = bre.tissue)
--) 
--select cs.*, case when bts.tissue is not null then bts.tissue else 'other' end as destination_tissue from contamined_samples cs
--left join bioqc_tissues_signatures bts on bts.signature = cs.signature
--order by tissue 
