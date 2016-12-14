------------------------------------------------------------
-- materialized view wit htissues relevant for the BioQC GEO study. 
-- could be replaced with a table where only relevant BioQC results
-- are imported in the first place. 
--
-- working on the full BioQC result table is infeasible for performance
-- reasons. 
-------------------------------------------------------------

create materialized view bioqc_res_tissue
parallel 16
build immediate
refresh force
on demand
as 
  with relevant_signatures as (
      select distinct bs.source from bioqc_tissue_set bts
      join bioqc_signatures bs on bs.id = bts.signature
  )
  select /*+ parallel(16)  */  br.* 
  from bioqc_res br
  join bioqc_signatures bs on bs.id = br.signature
  where bs.source in (
    --'gtex_ngs_0.85_5.gmt', 'exp.tissuemark.affy.roche.symbols.gmt'
    select  /*+ CARDINALITY(relevant_signatures, 2) */ source from relevant_signatures
  );
create /*+ parallel(16) */ index bioqc_res_tissue_gsm on bioqc_res_tissue(gsm);
create /*+ parallel(16) */ index bioqc_res_tissue_signature on bioqc_res_tissue(signature);
