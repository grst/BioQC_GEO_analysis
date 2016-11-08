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
create table bioqc_signatures(id serial primary key, 
    name varchar(255) not null,
    signature_set varchar(255) not null,
    description text null,
    unique(name, signature_set));

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

drop table if exists bioqc_genenames cascade;
create table bioqc_genenames(hgnc varchar(20) not null primary key,
    symbol varchar(80) not null unique,
    name text,
    status text,
    previous_symbols text,
    synonyms text,
    chromosome varchar(80),
    accession_numbers text,
    refseq_ids text);

drop table if exists bioqc_signatures_genes;
create table bioqc_signatures_genes(signature id int references bioqc_signatures(id),
    gene varchar(80) references bioqc_genenames(symbol),
    primary key(signature, gene));

