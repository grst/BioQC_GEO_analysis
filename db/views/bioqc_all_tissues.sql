create materialized view bioqc_all_tissues
build immediate
refresh force
on demand
as 
  with all_tissues as (
    select /*+ parallel(16) */  lower(tissue_orig) as tissue
                              , gsm as id
    from bioqc_gsm
    union 
    select /*+ parallel(16) */ lower(tissue_or_cell_type) as tissue
                             , cast(experiment_name as varchar(15)) as id
    from udis_meta
  )
  select /*+ parallel(16) */ tissue, count(id) as cnt
  from all_tissues
  group by tissue
  order by cnt desc;