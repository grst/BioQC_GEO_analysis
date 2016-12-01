/** Table Tissues: holds manually selected tissues **/
create table bioqc_tissues(id varchar2(80) not null primary key) 
  tablespace srslight_d;
  
/** map hetergenous tissue annotation from geo to unified tissue name **/
create table bioqc_normalize_tissues(tissue_orig varchar(1000) not null
                                   , tissue varchar(80) not null 
                                       references bioqc_tissues(id))
  tablespace srslight_d;

create table bioqc_signatures(id number(10) not null primary key
                            , name varchar2(255) not null          -- name of the signature in gmt
                            , source varchar2(255) not null        -- gmt filename
                            , description clob null                -- description in gmt
                            , gene_symbols clob null               -- comma separated list of gene symbols in gmt
                            , constraint uq_signatures
                                unique(name, source)
) tablespace srslight_d;

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
                                , primary key(tgroup, signature, tgroup, tissue_set)
    
) tablespace srslight_d;

/** for inserting tissue sets **/
create global temporary table bioqc_tmp_tissue_set (
    signature_name varchar2(255) not null 
  , signature_source varchar2(255) not null
  , tissue varchar(80) not null
  , tgroup varchar(80) not null
  , tissue_set varchar(80) not null
) on commit preserve rows 

create table bioqc_res(gsm varchar2(10) not null 
                        references bioqc_gsm(gsm)
                     , signature number(10) not null
                        references bioqc_signatures(id)
                        on delete cascade
                     , pvalue binary_double
                     , primary key(gsm, signature)
) tablespace srslight_d;

