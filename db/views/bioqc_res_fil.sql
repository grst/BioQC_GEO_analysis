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
as (
  select /*+ parallel(32) */  distinct  br.gsm
                                      , br.signature
                                      , br.pvalue 
                                      , bg.tissue
                                      , bg.organism_ch1 as organism
                                      , cast(cast(regexp_substr(submission_date, '^(\d{4})-.*', 1, 1, NULL, 1) as varchar2(4)) as NUMBER(4)) as year
                                      , cast(regexp_substr(contact, '.*Country:(.*?);.*?$', 1, 1, NULL, 1) as varchar2(100)) as country
  from bioqc_res br
  join bioqc_gsm bg on bg.gsm = br.gsm
  where bg.tissue is not null and bg.tissue != 'other'
  and channel_count = 1
  and organism_ch1 in ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')
);
create /*+ parallel(32) */ index gsm on bioqc_res_fil(gsm);
create /*+ parallel(32) */ index signature on bioqc_res_fil(signature);
