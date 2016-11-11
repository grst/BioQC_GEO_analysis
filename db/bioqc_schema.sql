/** Table Tissues: holds manually selected tissues **/
drop table if exists bioqc_tissues cascade;
create table bioqc_tissues(id varchar(80) not null primary key);

drop table if exists bioqc_signatures cascade;
create table bioqc_signatures(id varchar(80) not null primary key,
    description text null,
    gene_symbols text null);

drop table if exists bioqc_tissues_signatures cascade;
create table bioqc_tissues_signatures(tissue varchar(80) references bioqc_tissues(id),
    signature varchar(80) references bioqc_signatures(id),
    primary key(tissue, signature));

drop table if exists bioqc_res cascade;
create table bioqc_res(gse varchar(10) not null references gse(gse), 
    gsm varchar(10) not null references gsm(gsm),
    signature varchar(80) references bioqc_signatures(id),
    pvalue double precision,
    primary key(gse, gsm, signature));

