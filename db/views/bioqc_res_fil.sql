------------------------------------------------------------
-- filtered bioqc results that can be used for analysis
--   * tissue annotated (by running bioqc on GSE instead of the preselected
--     GSM we might introduce GSM without annotated tissue
--   * channel_count = 1 (remove obsolete studies, not considered in original
--     filtering
--   * Limit organisms to human, rat and mouse (most abundant ones) 
--
-- additionally, we add the parsed meta information year and country. 
--
-- remove gse annotation as some samples might be referenced in multiple
-- series and would therefore be duplicated. This is not desired. 
--
-- same samples in different studies indeed produce the same
-- bioqc result, which can be tested with the following query: 
-- 
-- ```
-- with duplicates as (
--   select /*+ parallel(16) */  gsm
--                             , count(distinct gse)
--   from bioqc_res
--   group by gsm 
--   having count(distinct gse) > 1
-- )
-- select /*+ parallel(16) */ *
-- from bioqc_res
-- join duplicates on duplicates.gsm = bioqc_res.gsm
-- order by signature, bioqc_res.gsm
-- ```
-------------------------------------------------------------

create materialized view bioqc_res_fil
parallel 16
build immediate
refresh force
on demand
as 
  with relevant_signatures as (
      select distinct bs.source from bioqc_tissue_set bts
      join bioqc_signatures bs on bs.id = bts.signature
  )
  select /*+ parallel(16)  */  distinct  br.gsm
                                      , br.signature
                                      , bs.name as signature_name
                                      , br.pvalue 
                                      , bnt.tissue
                                      , bg.organism_ch1 as organism
                                      , cast(
                                          cast(
                                            regexp_substr(submission_date, '^(\d{4})-.*', 1, 1, NULL, 1) 
                                            as varchar2(4)
                                          )
                                          as NUMBER(4)
                                        ) as year
                                      , cast(
                                          TRIM(BOTH from
                                               regexp_substr(contact, 'Country:(.*?)(;.*)?$', 1, 1, NULL, 1) 
                                          )
                                          as varchar2(100)
                                        ) as country
  from bioqc_res br
  join bioqc_gsm bg on bg.gsm = br.gsm
  join bioqc_normalize_tissues bnt on bnt.tissue_orig = lower(bg.tissue_orig)
  join bioqc_signatures bs on bs.id = br.signature
  where bs.source in (
   -- TODO: this is hardcoded and evil for performance reasons... need to correct this at a certain point. 
    'gtex_ngs_0.85_5.gmt', 'exp.tissuemark.affy.roche.symbols.gmt'
    --select  /*+ CARDINALITY(relevant_signatures, 2) */ source from relevant_signatures
  )
  and channel_count = 1
  and bg.organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus');
create /*+ parallel(16) */ index bioqc_res_fil_gsm on bioqc_res_fil(gsm);
create /*+ parallel(16) */ index bioqc_res_fil_signature on bioqc_res_fil(signature);
create index bioqc_res_fil_tissue on bioqc_res_fil(tissue);
