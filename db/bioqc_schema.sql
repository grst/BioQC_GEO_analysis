-- alter table bioqc_gsm add column tissue_orig
create index /*+ parallel(16) */ bioqc_gsm_tissue_orig
  on bioqc_gsm(lower(tissue_orig));

/** Table Tissues: holds manually selected tissues **/
create table bioqc_tissues(id varchar2(80) not null primary key) 
  tablespace srslight_d;
  
/** map hetergenous tissue annotation from geo to unified tissue name **/
create table bioqc_normalize_tissues(tissue_orig varchar(1000) not null primary key
                                   , tissue varchar(80) not null 
                                       references bioqc_tissues(id))
  tablespace srslight_d;
create index bioqc_normalize_tissues_tissue on bioqc_normalize_tissues(tissue);
  

create table bioqc_signatures(id number(10) not null primary key
                            , name varchar2(255) not null          -- name of the signature in gmt
                            , source varchar2(255) not null        -- gmt filename
                            , description clob null                -- description in gmt
                            , gene_symbols clob null               -- comma separated list of gene symbols in gmt
                            , constraint uq_signatures
                                unique(name, source)
) tablespace srslight_d;

create index bioqc_signatures_source on bioqc_signatures(source);

-- auto increment for bioqc_signatures
create sequence sig_seq start with 1;

create or replace trigger sig_bir
before insert on bioqc_signatures
for each row
begin 
  select sig_seq.NEXTVAL
  into :new.id
  from dual;
end;

create table bioqc_tissue_set( signature number(10) not null 
                                    references bioqc_signatures(id)
                                , tissue varchar2(80) not null
                                    references bioqc_tissues(id)
                                , tgroup varchar(80) not null
                                , tissue_set varchar(80) not null                  
                                , primary key(signature, tissue, tissue_set)
    
) tablespace srslight_d;
create index bioqc_tissue_set_tgroup on bioqc_tissue_set(tgroup); 
create index bioqc_tissue_set_tissue_set on bioqc_tissue_set(tissue_set); 
create index bioqc_tissue_set_tissue on bioqc_tissue_set(tissue); 
create index bioqc_tissue_set_signature on bioqc_tissue_set(signature); 
create index bioqc_tissue_set_tgts on bioqc_tissue_set(tgroup, tissue_set); 

/** for inserting tissue sets **/
create global temporary table bioqc_tmp_tissue_set (
    signature_name varchar2(255) not null 
  , signature_source varchar2(255) not null
  , tgroup varchar(80) not null
  , tissue varchar(80) not null
  , tissue_set varchar(80) not null
) on commit preserve rows 

CREATE global temporary TABLE bioqc_tmp_gse_gpl 
( gse varchar2(15),
	gpl varchar2(15),
  study_min float,
  study_25 float,
  study_median float, 
  study_mean float, 
  study_75 float,
  study_max float
) on commit preserve rows;

create table bioqc_res(gsm varchar2(10) not null 
                        references bioqc_gsm(gsm)
                     , signature number(10) not null
                        references bioqc_signatures(id)
                        on delete cascade
                     , pvalue binary_double
                     , primary key(gsm, signature)
) tablespace srslight_d;

-- contain all samples on which we successfully ran bioqc.
-- serves as background for contamination analysis
create table bioqc_bioqc_success(gsm varchar2(10) not null references bioqc_gsm(gsm) primary key)

