----------------------------------------------------------------------
-- ADD FOREIGN KEYS and FIX CONSTRAINTS
--
-- some tables are incomplete, we will add the ids with 
-- null in all other columns s.t. the constrainst are fulfilled. 
----------------------------------------------------------------------

-- works
alter table bioqc_gse_gpl add foreign key(gse) references bioqc_gse(gse);
alter table bioqc_gse_gpl add foreign key(gpl) references bioqc_gpl(gpl);

alter table bioqc_gse_gsm add foreign key(gse) references bioqc_gse(gse);

-- there are problems in:
-- gsm.gpl references gpl.gpl
-- gse_gsm.gsm references gsm.gsm
-- smatrix.gpl references gpl.gpl
-- smatrix.gse references gse.gse

-- fix GPL
insert into bioqc_gpl(gpl)
    select distinct gpl from bioqc_gsm where not exists(
        select * from bioqc_gpl where bioqc_gpl.gpl = bioqc_gsm.gpl);
insert into bioqc_gpl(gpl)
    select distinct gpl from bioqc_smatrix where not exists(
        select * from bioqc_gpl where bioqc_gpl.gpl = bioqc_smatrix.gpl); 

alter table bioqc_gsm add foreign key(gpl) references bioqc_gpl(gpl);
alter table bioqc_smatrix add foreign key(gpl) references bioqc_gpl(gpl);

-- fix GSM (6713 rows)  
insert into bioqc_gsm(gsm) 
    select distinct gsm from bioqc_gse_gsm where not exists( 
        select * from bioqc_gsm where bioqc_gsm.gsm = bioqc_gse_gsm.gsm); 

alter table bioqc_gse_gsm add foreign key(gsm) references bioqc_gsm(gsm); 

-- fix GSE  
insert into bioqc_gse(gse) 
    select distinct gse from bioqc_sMatrix where not exists( 
        select * from bioqc_gse where bioqc_gse.gse = bioqc_sMatrix.gse); 

alter table bioqc_smatrix add foreign key(gse) references bioqc_gse(gse); 
alter table bioqc_gds add foreign key(gpl) references bioqc_gpl(gpl);


------------------------------------------------
-- ADD ADDITIONAL COLUMNS
--
-- Add columns that will contain information 
-- we obtain from other sources
------------------------------------------------

alter table bioqc_gsm add (tissue varchar2(80) references bioqc_tissues("ID"));
alter table bioqc_gsm add (tissue_orig clob);

alter table bioqc_gpl add (has_annot number(1));
