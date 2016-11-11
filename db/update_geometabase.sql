-------------------------
-- ADD PRIMARY KEYS
-------------------------

alter table gds add primary key(gds);
alter table gds_subset add primary key("Name");
alter table geoconvert add primary key(from_acc, to_acc);
alter table geodb_column_desc add primary key("TableName", "FieldName");
alter table gpl add primary key(gpl);
alter table gse add primary key(gse);
alter table gsm add primary key(gsm);
alter table gse_gpl add primary key(gse, gpl);
alter table gse_gsm add primary key(gse, gsm);
alter table metainfo add primary key(name);
alter table smatrix add primary key("sMatrix");

----------------------------------------------------------------------
-- ADD FOREIGN KEYS and FIX CONSTRAINTS
--
-- some tables are incomplete, we will add the ids with 
-- null in all other columns s.t. the constrainst are fulfilled. 
----------------------------------------------------------------------

-- works
alter table gse_gpl add foreign key(gse) references gse(gse);
alter table gse_gpl add foreign key(gpl) references gpl(gpl);

alter table gse_gsm add foreign key(gse) references gse(gse);

-- there are problems in:
-- gsm.gpl references gpl.gpl
-- gse_gsm.gsm references gsm.gsm
-- smatrix.gpl references gpl.gpl
-- smatrix.gse references gse.gse

-- fix GPL
insert into gpl(gpl)
    select distinct gpl from gsm where not exists(
        select * from gpl where gpl.gpl = gsm.gpl);
insert into gpl(gpl)
    select distinct gpl from sMatrix where not exists(
        select * from gpl where gpl.gpl = sMatrix.gpl); 

alter table gsm add foreign key(gpl) references gpl(gpl);
alter table smatrix add foreign key(gpl) references gpl(gpl);

-- fix GSM (6713 rows)  
insert into gsm(gsm) 
    select distinct gsm from gse_gsm where not exists( 
        select * from gsm where gsm.gsm = gse_gsm.gsm); 

alter table gse_gsm add foreign key(gsm) references gsm(gsm); 

-- fix GSE  
insert into gse(gse) 
    select distinct gse from sMatrix where not exists( 
        select * from gse where gse.gse = sMatrix.gse); 

alter table smatrix add foreign key(gse) references gse(gse); 
alter table gds add foreign key(gpl) references gpl(gpl);


------------------------------------------------
-- ADD ADDITIONAL COLUMNS
--
-- Add columns that will contain information 
-- we obtain from other sources
------------------------------------------------

alter table gsm add column tissue varchar(80) references bioqc_tissues(id);
alter table gsm add column tissue_orig text;

alter table gpl add column has_annot boolean;
