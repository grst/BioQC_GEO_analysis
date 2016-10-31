/** Table Tissues: holds manually selected tissues **/
drop table if exists bioqc_tissues cascade;
create table bioqc_tissues(id varchar(80) not null primary key);

/** 
 * Table GEO Samples: holds additional meta information that is not
 * contained in GEOmetabase (e.g. the processed tissue information)
 */
drop table if exists bioqc_gsm cascade;
create table bioqc_gsm(gsm varchar(10) primary key references gsm(gsm),
    tissue varchar(80) references bioqc_tissues(id),
    tissue_orig text);

drop table if exists bioqc_signatures cascade;
create table bioqc_signatures(id varchar(80) not null primary key,
    description text null,
    gene_symbols text null);

drop table if exists bioqc_tissues_signatures cascade;
create table bioqc_tissues_signatures(tissue varchar(80) references bioqc_tissues(id),
    signature varchar(80) references bioqc_signatures(id),
    primary key(tissue, signature));

drop table if exists bioqc_res cascade;
create table bioqc_res(gsm varchar(10) not null references gsm(gsm),
    signature varchar(80) references bioqc_signatures(id),
    pvalue double precision,
    primary key(gsm, signature));

