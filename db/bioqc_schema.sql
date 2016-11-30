/** Table Tissues: holds manually selected tissues **/
create table bioqc_tissues(id varchar2(80) not null primary key) 
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

create table bioqc_signature_set( signature number(10) not null 
                                    references bioqc_signatures(id)
                                , tissue varchar2(80) not null
                                    references bioqc_tissues(id)
                                , tgroup varchar(80) not null
                                , signature_set varchar(80) not null                  
                                , constraint pk_bioqc_signature_set
                                    primary key(tissue, signature, tgroup)
    
) tablespace srslight_d;

create table bioqc_res(gse varchar2(10) not null 
                     , gsm varchar2(10) not null
                     , signature number(10) not null
                     , pvalue binary_double
                     , primary key(gse, gsm, signature)
                     , constraint fk_gse
                        foreign key (gse)
                        references bioqc_gse(gse)
                     , constraint fk_gsm
                        foreign key (gsm)
                        references bioqc_gsm(gsm)
                     , constraint fk_bioqc_res_signature
                        foreign key (signature)
                        references bioqc_signatures(id)
) tablespace srslight_d;

