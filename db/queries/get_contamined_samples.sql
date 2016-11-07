-- using group by having count, works on filtered table, but not as flexible. 
select gse, gsm, tissue, signature, min(enrichment_ratio) from bioqc_res_ext2 bre
where tissue = 'adipose'
group by gse, gsm, tissue, signature
having count(exp_sig) = (select count(*) from bioqc_tissues_signatures where tissue = bre.tissue)

