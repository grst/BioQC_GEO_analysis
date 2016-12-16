# Tissue Migration




Now, detect samples we defined as contamined:
* highly siginificant BioQC scores (min score per sample >6)
* signature outside 'expected' signature that is at least 6 higher than the highest score in expected signatures. 
* (outlier in study)

We can find contamined studies with SQL only. 
To speed up things, we create a temporary table (rather a materialized view, but that's 
not available on our old postgres server) filtering for samples, where a non-expected
signature exceeds all expected signatures by a certain threshold (BioQC score > 6, aka. 'enrichment-ratio'). This table is created in `db/queries/mk-extended_bioqc_res_table2.sql`. 

We create a view `contamined_samples` that contains all contamined samples in `db/queries/get_contamined_samples.sql`. 

Second, we retrieve the information we need for the migration chart from that temporary table. 

```r
query = "
select origin
     , destination
     , count(gsm) as \"count\"
from bioqc_tissue_migration 
where enrichment_ratio > 4
  and rk = 1
  and tissue_set = 'gtex_solid'
group by origin, destination
order by origin, destination"
migration = dbGetQuery(mydb, query)
set.seed(42)
col = rand_color(length(unique(migration$ORIGIN)))
chordDiagram(migration, grid.col=col, annotationTrack="grid", preAllocateTracks=1)
circos.trackPlotRegion(track.index=1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
```


```r
chordDiagram(migration[migration$ORIGIN != 'blood' & migration$DESTINATION!='blood',])
```

Contamination in total: 


We can also compose statistics about the origin of samples

By time

```r
byYear = "
with contam_stats as (
  select bcs.*
       , case when enrichment_ratio > 4
           then 1
           else null
         end as is_contam 
  from bioqc_contam_stats bcs
)
select  year
       , tgroup
       , count(gsm) as total
       , count(is_contam) as contamined
       , count(is_contam)/cast(count(gsm) as float) as ratio
from contam_stats
where tissue_set = 'gtex_solid'
group by year, tgroup
order by ratio desc
"

year = dbGetQuery(mydb, byYear)

ggplot(data=year, aes(x=YEAR, y=RATIO)) + stat_summary(fun.y="mean", geom="bar")
```


By country
(that's nice figures, but it's not very representative... tissue-bias, statistical significance, ...)

```r
byCountry = "
with contam_stats as (
  select bcs.*
       , case when enrichment_ratio > 4
           then 1
           else null
         end as is_contam 
  from bioqc_contam_stats bcs
)
select  country
       , tgroup
       , count(gsm) as total
       , count(is_contam) as contamined
       , count(is_contam)/cast(count(gsm) as float) as ratio
from contam_stats
where tissue_set = 'gtex_solid'
group by country, tgroup
order by ratio desc
"

country = dbGetQuery(mydb, byCountry)
country.2 = country[country$TGROUP == 'liver',]
#country.2 = ddply(country[country$total > 2000,], ~country,summarise,ratio=mean(ratio))

spdf <- joinCountryData2Map(country.2, joinCode="NAME", nameJoinColumn="COUNTRY")
mapCountryData(spdf, nameColumnToPlot="RATIO", catMethod="fixedWidth")
```

By organism

```r
byOrganism = "
with contam_stats as (
  select bcs.*
       , case when enrichment_ratio > 4
           then 1
           else null
         end as is_contam 
  from bioqc_contam_stats bcs
)
select  organism
       , tgroup
       , count(gsm) as total
       , count(is_contam) as contamined
       , count(is_contam)/cast(count(gsm) as float) as ratio
from contam_stats
where tissue_set = 'gtex_solid'
group by organism, tgroup
order by ratio desc
"

organism = dbGetQuery(mydb, byOrganism)
tissues = c("blood", "brain", "kidney", "liver")
organism.2 = organism[organism$TGROUP %in% tissues, ]

#ggplot(data=organism.2, aes(x=ORGANISM, y=RATIO)) + stat_summary(fun.y="mean", geom="bar")
ggplot(data=organism.2, aes(x=ORGANISM, y=RATIO)) + 
  geom_boxplot() +
  geom_point(aes(color=TGROUP))
```

