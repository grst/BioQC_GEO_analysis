# Results [PRELIMINARY]




```r
sql = "
select /*+ parallel(16) */ *
from bioqc_contamination
where tissue_set = 'gtex_solid'
"
res = data.table(dbGetQuery(mydb, sql))
res = res[,MIN_FOUND_PVALUE:=as.numeric(MIN_FOUND_PVALUE)]
res = res[,MIN_EXP_PVALUE:=as.numeric(MIN_EXP_PVALUE)]
res = res[,EXP_SCORE:=absLog10p(MIN_EXP_PVALUE)]
res = res[,FOUND_SCORE:=absLog10p(MIN_FOUND_PVALUE)]
res = res[,total_tissue:=length(unique(GSM)), by=c("TGROUP")]
#res[,cnt_per_group:=length(unique(GSM)), by=c("GPL", "TGROUP")]
#res = res[cnt_per_group >= 100]
pbonf = 0.05 / nrow(res)
```


## Detection of tissues
In general, a signature finds its corresponding tissue: 


```r
ggplot(res, aes(x=FOUND_TGROUP, y=FOUND_SCORE)) + 
  geom_boxplot() + facet_wrap(~TGROUP) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-2-1.png" width="1152" style="display:block; margin: auto" style="display: block; margin: auto;" />

Although, for each tissue there are samples in which the signature is not enriched (pvalue >= 0.05) at a single sample level: 


```r
matching = res[TGROUP == FOUND_TGROUP]
not_enriched = data.table(sqldf("select tgroup
                            , sum(case when MIN_FOUND_PVALUE >= 0.05 then 1 else 0 end) as not_enriched
                            , count(distinct gsm) as total 
                      from matching
                      group by tgroup"))
```

```
## Loading required package: tcltk
```

```
## Warning: Quoted identifiers should have class SQL, use DBI::SQL() if the
## caller performs the quoting.
```

```r
not_enriched[, ratio:=not_enriched/total]
```

```
##             TGROUP not_enriched total      ratio
## 1:           blood         1331 19393 0.06863301
## 2:           brain          504  7624 0.06610703
## 3:           heart           36  2402 0.01498751
## 4:          kidney          345  3540 0.09745763
## 5:           liver          315 14019 0.02246951
## 6:        pancreas          297   586 0.50682594
## 7: skeletal muscle           98  1830 0.05355191
## 8:            skin          321  1927 0.16658018
## 9:          testis          418   627 0.66666667
```

```r
ggplot(not_enriched, aes(x=TGROUP, y=ratio)) + geom_bar(stat="identity")
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display:block; margin: auto" style="display: block; margin: auto;" />


And, there are outliers (=tissue heterogeniety). 

## Detection of outliers

### Sample-against-sample method

1) it is not a valid approach to compare samples against other samples. Even if stratified by platform, different experimental conditions and data processing steps could influence the BioQC score in a way that it would fall outside a confidence interval.
2) taking (1) into account we could detect outliers at a study level. This can be trivially done by computing a z-score on the BioQC scores. The statistically correct approach to check if a signature is *more enriched* in one sample than another would be as follows: 

For two samples $S_1$ and $S_2$ and a signature $K$:
$$
H_0: \text{genes in K are not more enriched in } S_1 \text{ than in } S_2 \\
H_1: \text{genes in K are more enriched in } S_1 \text{ than in } S_2 
$$ 
which translates in the context of the WMW test into
$$
H_1: \text{The median rank of signature genes in } S_1 \text{ is higher than the median rank of signature genes } S_2 
$$

Let $U_1$ and $U_2$ be the WMW test statistics for $S_1$ and $S_2$ respectively, with 

$$
U_1 \sim \mathcal{N}(\mu_1, \sigma_1^2) \\
U_2 \sim \mathcal{N}(\mu_2, \sigma_2^2)
$$

We define the random variable 
$$
\Delta = U_1 - U_2 
$$

Which is distributed, assuming independence of the two samples, according to
$$
\Delta \sim \mathcal{N}(\mu_1 - \mu_2, \sigma_1^2 - \sigma_2^2).
$$

So that we now can easily test our hypothesis: 
$$
H_1: \Delta > 0. 
$$

### Within-Sample Method

Another approach is to check, which signature scores highest in a sample. Like that we classify a sample based on the signature, independent of the original annotation. As shown in the [signature validation](#Validating Tissue Signatures) this works reliably with a sensitivity and specificity of 1.0 for the selected subset of signatures. 
Acknowledging that we have tested the signatures only on a small set of samples, we can add an additional safety margin, by requiring that the p-value exceeds the p-value of the 'expected signature' at least by 2 orders of magnitude. 

The downside of this approach is, that it will only detect heavily contaminated or mislabeled samples. Also, some signatures have a higher score than others in general, which further limits our detection limit. For instance, it will be relatively easy to detect liver contamination in pancreas, while it is hard to detect pancreas contamination in liver. 

We filter the results according to this criterion: 

```r
res = res[,score_ratio:=FOUND_SCORE - EXP_SCORE]
fil = res[TGROUP != FOUND_TGROUP & score_ratio > 2]

#### exclude 'double contaminations'
fil = fil[, rk:=frankv(FOUND_SCORE, order=-1), by=c("GSM", "GPL", "TGROUP")]
fil = fil[rk==1]
```


And display the 'contamination matrix' 

```r
aggr = sqldf('select TGROUP, FOUND_TGROUP, count(GSM) as cnt from fil group by TGROUP, FOUND_TGROUP')
ggplot(aggr, aes(y=TGROUP, x=FOUND_TGROUP)) +
  geom_tile(aes(fill=cnt)) +
  geom_text(aes(label=cnt))
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display:block; margin: auto" style="display: block; margin: auto;" />

The same can be visualized as 'migration chart'

```r
set.seed(42)
col = rand_color(length(unique(aggr$TGROUP)))
names(col) = unique(aggr$TGROUP)
chordDiagram(aggr, grid.col=col, annotationTrack="grid", preAllocateTracks=1)
circos.trackPlotRegion(track.index=1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=2)
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-6-1.png" width="480" style="display:block; margin: auto" style="display: block; margin: auto;" />

The fraction of contaminated samples per tissue: 

```r
per_tissue = sqldf("select tgroup
                  , total_tissue as TOTAL
                   , count(gsm) as CONTAMINATED
                   , count(gsm)/cast(total_tissue as float) as RATIO
                   from fil group by tgroup")
# png(filename = "contamination.png", width = 1200, height=600)
ggplot(data=per_tissue, aes(x=TGROUP, y=RATIO)) +
  geom_bar(stat="identity") + 
  xlab("tissue") + 
  ylab("fraction contaminated")
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display:block; margin: auto" style="display: block; margin: auto;" />

These are actually so few examples, that we can have a look at them manually: 

```r
sql = "
select * from bioqc_gse_gsm
"
all_gse = data.table(dbGetQuery(mydb, sql))
setkey(all_gse, GSM)
contam_studies = merge(fil, all_gse, all=FALSE, by=c("GSM"))
datatable(sqldf("select GSE, GSM, GPL, TGROUP, FOUND_TGROUP, ORGANISM, EXP_SCORE, FOUND_SCORE, score_ratio 
       from contam_studies
       order by gse asc, score_ratio desc"))
```

<!--html_preserve--><div id="htmlwidget-158363e772d6bccf63af" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-158363e772d6bccf63af">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276"],["GSE11990","GSE11990","GSE14668","GSE14668","GSE14668","GSE14668","GSE16595","GSE17059","GSE17060","GSE17183","GSE19044","GSE19044","GSE19044","GSE20524","GSE20524","GSE20524","GSE20524","GSE20524","GSE20524","GSE20524","GSE21691","GSE22417","GSE22417","GSE22417","GSE22417","GSE24637","GSE25630","GSE25632","GSE25846","GSE25846","GSE25846","GSE26600","GSE27045","GSE27045","GSE31973","GSE31973","GSE31973","GSE31973","GSE31973","GSE31973","GSE31973","GSE31973","GSE32993","GSE36443","GSE36443","GSE36443","GSE36452","GSE36452","GSE36452","GSE36618","GSE36618","GSE36618","GSE39759","GSE39759","GSE39759","GSE39759","GSE39759","GSE39759","GSE39759","GSE39759","GSE40225","GSE40225","GSE40435","GSE40435","GSE43339","GSE45513","GSE45551","GSE45551","GSE45551","GSE45551","GSE45551","GSE45551","GSE45551","GSE45551","GSE45551","GSE49000","GSE50827","GSE50949","GSE50949","GSE50949","GSE50949","GSE50949","GSE50949","GSE52004","GSE52004","GSE52004","GSE52141","GSE52141","GSE52141","GSE52141","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE52403","GSE53605","GSE53605","GSE53605","GSE53605","GSE53605","GSE53605","GSE54015","GSE54480","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE54674","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE55998","GSE60506","GSE60506","GSE62029","GSE62029","GSE62029","GSE62037","GSE62037","GSE62037","GSE62732","GSE62732","GSE62732","GSE63362","GSE63362","GSE63362","GSE63362","GSE63362","GSE63362","GSE63362","GSE63793","GSE66420","GSE66420","GSE66420","GSE66420","GSE66420","GSE67225","GSE67225","GSE67225","GSE67225","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67278","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE67285","GSE72925","GSE73759","GSE73759","GSE75225","GSE75225","GSE80427","GSE80431","GSE83562","GSE83562","GSE83562"],["GSM303382","GSM303379","GSM366059","GSM366060","GSM366057","GSM366058","GSM417047","GSM426799","GSM426799","GSM429412","GSM471358","GSM471355","GSM471354","GSM515585","GSM515563","GSM515592","GSM515593","GSM515581","GSM515565","GSM515588","GSM541266","GSM570065","GSM570062","GSM570063","GSM570064","GSM607586","GSM629753","GSM629753","GSM634838","GSM634840","GSM634842","GSM654820","GSM667541","GSM667544","GSM791842","GSM791841","GSM791844","GSM791843","GSM791838","GSM791839","GSM791837","GSM791840","GSM817263","GSM893619","GSM893620","GSM893621","GSM893619","GSM893620","GSM893621","GSM897564","GSM897563","GSM897562","GSM978799","GSM978800","GSM978798","GSM978801","GSM978803","GSM978802","GSM978804","GSM978805","GSM988591","GSM988589","GSM993984","GSM994103","GSM1060691","GSM1105884","GSM1109140","GSM1109139","GSM1109130","GSM1109132","GSM1109125","GSM1109137","GSM1109144","GSM1109128","GSM1109129","GSM1191961","GSM1230445","GSM1233153","GSM1233151","GSM1233152","GSM1233149","GSM1233148","GSM1233150","GSM1257094","GSM1257092","GSM1257093","GSM1260251","GSM1260250","GSM1260252","GSM1260253","GSM1265286","GSM1265051","GSM1265015","GSM1265042","GSM1265043","GSM1265063","GSM1265139","GSM1265090","GSM1265097","GSM1265050","GSM1265197","GSM1265289","GSM1265171","GSM1265162","GSM1265463","GSM1265149","GSM1265041","GSM1265163","GSM1265216","GSM1265030","GSM1265148","GSM1265099","GSM1265147","GSM1265374","GSM1265210","GSM1265470","GSM1265373","GSM1265140","GSM1265187","GSM1265082","GSM1265284","GSM1265195","GSM1265209","GSM1265040","GSM1265014","GSM1265332","GSM1265091","GSM1265161","GSM1265049","GSM1265146","GSM1265118","GSM1265196","GSM1265112","GSM1265244","GSM1265138","GSM1265189","GSM1265130","GSM1265070","GSM1296703","GSM1296684","GSM1296669","GSM1296667","GSM1296673","GSM1296668","GSM1305882","GSM1316482","GSM1321588","GSM1321615","GSM1321574","GSM1321613","GSM1321607","GSM1321625","GSM1321573","GSM1321610","GSM1321596","GSM1321611","GSM1321612","GSM1321609","GSM1321572","GSM1321618","GSM1321587","GSM1321605","GSM1321581","GSM1321575","GSM1321577","GSM1321599","GSM1321619","GSM1321622","GSM1321583","GSM1321586","GSM1321606","GSM1321582","GSM1321601","GSM1321608","GSM1321624","GSM1321584","GSM1321602","GSM1321585","GSM1321594","GSM1321617","GSM1321580","GSM1321576","GSM1321620","GSM1321616","GSM1321592","GSM1321597","GSM1321578","GSM1321593","GSM1321614","GSM1321598","GSM1321590","GSM1321604","GSM1321591","GSM1321579","GSM1350066","GSM1350054","GSM1350064","GSM1350063","GSM1350062","GSM1350069","GSM1350057","GSM1350065","GSM1350058","GSM1350061","GSM1350055","GSM1350056","GSM1350052","GSM1350068","GSM1350053","GSM1350060","GSM1350059","GSM1350067","GSM1481346","GSM1481345","GSM1518701","GSM1518699","GSM1518700","GSM1518701","GSM1518699","GSM1518700","GSM1532547","GSM1532546","GSM1532548","GSM1547455","GSM1547451","GSM1547458","GSM1547456","GSM1547450","GSM1547457","GSM1547452","GSM1557525","GSM1622284","GSM1622283","GSM1622287","GSM1622285","GSM1622281","GSM1642679","GSM1642678","GSM1642674","GSM1642675","GSM1643644","GSM1643645","GSM1643659","GSM1643646","GSM1643661","GSM1643660","GSM1643642","GSM1643640","GSM1643654","GSM1643643","GSM1643655","GSM1643658","GSM1643641","GSM1643644","GSM1643645","GSM1643659","GSM1643646","GSM1643661","GSM1643660","GSM1643642","GSM1643640","GSM1643654","GSM1643643","GSM1643655","GSM1643658","GSM1643641","GSM1874828","GSM1902266","GSM1902265","GSM1946372","GSM1946360","GSM2127137","GSM2127137","GSM2209654","GSM2209659","GSM2209657"],["GPL1261","GPL1261","GPL570","GPL570","GPL570","GPL570","GPL1261","GPL6101","GPL6101","GPL570","GPL6887","GPL6887","GPL6887","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL6246","GPL8321","GPL8321","GPL8321","GPL8321","GPL1261","GPL6884","GPL6884","GPL6246","GPL6246","GPL6246","GPL1261","GPL7546","GPL7546","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL6480","GPL4134","GPL4134","GPL4134","GPL4134","GPL4134","GPL4134","GPL1261","GPL1261","GPL1261","GPL7202","GPL7202","GPL7202","GPL7202","GPL7202","GPL7202","GPL7202","GPL7202","GPL6246","GPL6246","GPL10558","GPL10558","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL6885","GPL10558","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL8321","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL11533","GPL10558","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL6885","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL571","GPL1261","GPL1261","GPL570","GPL570","GPL570","GPL570","GPL570","GPL570","GPL1261","GPL1261","GPL1261","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1355","GPL1261","GPL6246","GPL6246","GPL6246","GPL6246","GPL6246","GPL6887","GPL6887","GPL6887","GPL6887","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL1261","GPL570","GPL7202","GPL7202","GPL6246","GPL6246","GPL6246","GPL6246","GPL6887","GPL6887","GPL6887"],["skin","skin","liver","liver","liver","liver","skeletal muscle","liver","liver","liver","liver","liver","liver","blood","blood","blood","blood","blood","blood","blood","skin","skin","skin","skin","skin","liver","brain","brain","brain","brain","brain","brain","kidney","kidney","liver","liver","liver","liver","liver","liver","liver","liver","skin","liver","liver","liver","liver","liver","liver","liver","liver","liver","brain","brain","brain","brain","brain","brain","brain","brain","pancreas","pancreas","kidney","kidney","liver","skin","skin","skin","skin","skin","skin","skin","skin","skin","skin","skeletal muscle","pancreas","liver","liver","liver","liver","liver","liver","kidney","kidney","kidney","liver","liver","liver","liver","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","kidney","kidney","kidney","kidney","kidney","kidney","kidney","skin","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","brain","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","heart","heart","liver","liver","liver","liver","liver","liver","kidney","kidney","kidney","skin","skin","skin","skin","skin","skin","skin","liver","brain","brain","brain","brain","brain","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","kidney","skin","skin","liver","liver","skin","skin","brain","brain","brain"],["skeletal muscle","skeletal muscle","blood","blood","blood","blood","liver","kidney","kidney","blood","blood","blood","blood","liver","liver","liver","liver","liver","liver","liver","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","blood","blood","blood","blood","blood","liver","blood","blood","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","blood","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","liver","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","liver","skeletal muscle","skeletal muscle","blood","blood","blood","blood","liver","blood","liver","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","kidney","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","skeletal muscle","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","skeletal muscle","skeletal muscle","skeletal muscle","blood","blood","skeletal muscle","skeletal muscle","blood","blood","blood"],["Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus","Rattus norvegicus","Rattus norvegicus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus","Mus musculus","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus"],[12.666296788985,11.559748458809,9.83002085410823e-06,1.05806337630452e-09,4.67255230646243e-12,9.95107835072207e-11,2.24320413597406,6.42798315700797,6.42798315700797,0,0.0434820918673718,2.27638264544473,2.72569810202282,5.32688296764215,11.8671833212037,9.54323207178514,7.88637274674847,8.32680458853306,14.4714548347203,19.7354927988327,23.2192142663006,0.83152650844089,2.76295764317846,3.3604960486783,2.75270932085954,1.08015573167773e-10,0.000289696446560811,0.000289696446560811,0,0,0,0,0.151751028091325,0.309259278744973,1.22986473263686e-06,1.23583924581623e-06,1.23367285847568e-06,1.2284267565254e-06,1.24578918754133e-06,1.24516507634207e-06,1.23480717897948e-06,1.23552953911506e-06,3.47157887959185e-15,2.95151700282593,3.05087956050962,9.86907168399409e-05,2.95151700282593,3.05087956050962,9.86907168399409e-05,0,0,0,0,0,0,0,0,0,3.03763151964287e-15,4.33947359948979e-16,2.03891710043836e-10,7.53477265992064e-13,0.00552105425594912,0.0260893759293382,0,7.1850489171737,9.26059739896217,9.36755594600998,16.1208124311808,13.26519832517,20.1093236255681,13.1125425693277,14.1126092853557,16.4733632193343,13.6410490175319,27.6837977917044,0.00363467727581904,0,0,0,0,0,0,0.227380572865747,0.381252876442088,0.265460281686316,0,0,0,0,5.85737522736227,6.85632158814127,7.05644186198086,7.78587319827022,7.13708552471862,10.0359413345282,8.81983868726053,8.96449181640129,12.2331406799886,8.5559745145947,4.17475732028469,9.30255824340348,13.8598952519207,14.2852885626299,12.8644231032967,9.6508860675157,11.1865422458743,10.7252371810493,11.0347194159204,11.9033044791993,9.43173574701637,11.6846618003792,12.7007962107156,13.084058843218,16.4092113629279,13.2972307803066,15.4982999315105,8.68486679563263,12.3496853030033,11.3985682804932,13.7538644610557,19.1232623867463,11.9389998904736,13.3840851882001,15.0327251389433,15.2393221743648,11.1657606569602,11.2710106049854,11.1079126683388,19.1788341276903,13.4626483314036,12.6541018478823,19.2528171202641,12.4777137072246,11.9978108063444,9.09733114626106,12.5739726910072,14.4816227802385,6.99205447116752e-10,1.99123978573153e-11,0.409551853841986,0.000381313179918665,0.0806576173442244,0.122158109296569,28.5326851458081,0.000103338069834087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.01825555158472e-10,0,2.84680987951163e-09,7.9995871851075e-10,0,1.98892539976661e-13,0,0,8.67894719897959e-16,7.7348706092308e-13,3.47157887959185e-15,4.02995861402943e-08,7.27725726772521e-09,1.16011824878323e-09,0,1.16757442736207e-10,0,0,0,4.88246947495362e-05,1.56128864978109e-09,5.10578317407993e-10,4.88246947495362e-05,1.56128864978109e-09,5.10578317407993e-10,0.000102290637942156,0.00139692381672708,1.4819345739522e-10,3.61167622994852,4.66265547402751,7.18289707980512,6.44963447647952,5.52854725776566,5.57790182338508,7.67751632562153,2.6806526333963,0,0,0,0,0,1.42436422417971e-10,5.14405398354468e-09,0.000304110698193556,0.000749643527812753,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18.0690075315515,18.3971635820766,0.0985386795320964,1.20438714609376e-11,24.1633056194392,24.1633056194392,0,0,0],[17.3262926720887,15.1200552721112,4.86012825069262,3.59408364601632,3.55423797187294,2.7271020024827,15.0712043821351,20.9322183829808,20.9322183829808,23.0546330592159,14.5333081523597,6.86474730433178,5.85894451285312,52.896041830597,44.6399381600343,38.0502457442144,36.3385238805086,34.3863718598295,29.7770251571763,21.7956637686174,28.512221948289,5.93895535162827,5.94101349000819,6.52104176306547,5.7924257107307,8.86132809344627,2.04558881703237,2.04558881703237,4.03834913966401,3.78396100979512,2.87589879450391,57.187564863265,2.38025006243373,2.52312855701228,15.8371245896006,15.8367377207018,15.8351129252567,15.8350355562561,15.8249019201219,15.8240511583956,15.8228910672898,15.8210350139705,10.7844733488977,51.539170054769,51.3442868560219,44.9315431827628,51.539170054769,51.3442868560219,44.9315431827628,19.1057598489041,13.6560249528921,13.1341701780478,10.561830976763,9.63739613780582,9.09041491177564,8.17930030929565,5.58877747151777,4.71399814914496,4.33417873177785,3.86629340025828,8.72394349111785,8.61140526727815,2.65648071963889,2.17094471430853,3.11640699050514,32.794631546185,19.4796934320364,18.090374321321,23.3717128947081,20.0837526715134,25.1642755296271,17.7200251060117,17.2841307712188,19.6174618930308,16.6306517368626,32.6856909134185,2.43218089817186,30.7190780280846,30.1207792736232,28.4962161488191,28.1439002010615,27.8137262864572,27.8113718397785,7.36400353566855,7.07456081955256,6.30297360049969,8.81951613865488,6.51521151539283,2.78973693455129,2.62745564169631,51.9632071797022,50.0170450777575,48.1982058052661,48.4020914224123,46.9160711105841,46.1885252714895,41.7059418442287,40.7300403161398,39.8812971070281,33.6757399496309,29.1555332001438,32.8648041192453,36.2516350383853,36.1853104543059,34.6291084015593,31.3213650285816,32.0801335185106,31.5067189760893,31.3428843946545,32.0400919957229,29.2796306733898,31.1237615311911,31.6370665064998,31.4419617969377,34.5798675934724,29.1517756064816,31.1771356587705,23.3028805126707,26.6807915893967,24.9330399932663,26.6754290565778,30.9921115432061,23.0968644063041,24.2294287711731,25.2093294737117,24.5873914548834,19.6887554835782,19.7638981603083,19.1664680985829,27.1389179767471,21.1983251798477,17.909448974912,24.06419560918,17.1174221585683,16.2553172498035,12.5992510483443,15.5033315801198,17.169053183942,21.17296973231,9.42328666353569,4.11034586817046,3.44578778124797,2.43291629249752,2.38192518993612,34.2379019554778,3.20437826709499,60.0362308806446,18.6015942045366,17.8639342436323,17.8434073714439,17.8126760727161,17.7950558141117,17.4545924703922,17.3071465749948,17.1936703619909,17.1681847748447,17.1186819642422,17.1116363453666,17.1102274030416,17.0332503959747,16.9763730444757,16.9537900972103,16.9171437492488,16.9019508067597,16.8623566324271,16.8217635160776,16.8175045483819,16.8050799566303,16.7695079288439,16.6942391737494,16.6564973735799,16.6502461322852,16.5656278687299,16.513720103301,16.5080500099678,16.4899401811273,16.3987198181955,16.3508007043865,16.2713040538615,16.1493281434008,16.1422202730942,16.0660502124752,16.0573258902284,15.987892400769,15.9352416438171,15.9176663288293,15.7064076576521,15.5468375298472,15.5466365062712,15.5359842016412,15.4400849516219,15.4339430193406,14.9876530007503,14.9636693680233,22.8567432425572,22.8040258782399,22.7752047640768,22.6007531033818,21.26084718456,21.0477373761922,20.9413832884132,20.2692149084912,19.9698646934514,19.8739349560792,19.6446304322621,19.4820400187451,18.6986827492858,17.6281539144207,16.0812073785277,15.9435573712599,13.9857104932406,9.31161927723081,4.06011169260699,3.08725473489223,4.8285629488342,3.62338475133462,2.71829234495848,4.8285629488342,3.62338475133462,2.71829234495848,3.89311048294548,3.38309693027673,2.87381372556116,14.2487218914694,13.8321282857128,15.3447402786978,14.5479129269877,13.4935994720086,11.2950012215201,10.9993226260937,19.4576393735276,3.7638280318425,3.51748478321809,3.46056427117438,2.50663986351731,2.33359184059253,22.9713556446406,22.7030741978745,8.96458147096435,8.70627015950135,12.5924794344816,11.9814215696202,10.8698804675621,10.1771522960867,7.93200317859503,7.81331144540585,4.71428189158595,4.14078220474963,3.69854391885719,3.48781329054908,2.44502954491056,2.29601777557942,2.240925083768,12.5924794344816,11.9814215696202,10.8698804675621,10.1771522960867,7.93200317859503,7.81331144540585,4.71428189158595,4.14078220474963,3.69854391885719,3.48781329054908,2.44502954491056,2.29601777557942,2.240925083768,30.8732891020125,24.6612950167218,22.9251245222056,2.78589591553828,2.29274023214574,26.5172193021324,26.5172193021324,2.22204068529081,2.21216723311021,2.16217287414067],[4.65999588310362,3.56030681330219,4.86011842067177,3.59408364495826,3.55423797186827,2.72710200238319,12.828000246161,14.5042352259728,14.5042352259728,23.0546330592159,14.4898260604924,4.58836465888705,3.13324641083031,47.5691588629549,32.7727548388306,28.5070136724292,28.4521511337601,26.0595672712965,15.305570322456,2.06017096978463,5.29300768198842,5.10742884318738,3.17805584682973,3.16054571438717,3.03971638987116,8.86132809333825,2.04529912058581,2.04529912058581,4.03834913966401,3.78396100979512,2.87589879450391,57.187564863265,2.22849903434241,2.21386927826731,15.8371233597359,15.8367364848625,15.8351116915839,15.8350343278293,15.8249006743327,15.8240499132305,15.8228898324826,15.821033778441,10.7844733488977,48.5876530519431,48.2934072955123,44.931444492046,48.5876530519431,48.2934072955123,44.931444492046,19.1057598489041,13.6560249528921,13.1341701780478,10.561830976763,9.63739613780582,9.09041491177564,8.17930030929565,5.58877747151777,4.71399814914496,4.33417873177785,3.86629340025828,8.72394349091396,8.6114052672774,2.65095966538294,2.14485533837919,3.11640699050514,25.6095826290113,10.2190960330742,8.72281837531099,7.25090046352738,6.81855434634334,5.0549519040589,4.60748253668395,3.17152148586307,3.14409867369654,2.98960271933071,5.00189312171408,2.42854622089605,30.7190780280846,30.1207792736232,28.4962161488191,28.1439002010615,27.8137262864572,27.8113718397785,7.13662296280281,6.69330794311047,6.03751331881337,8.81951613865488,6.51521151539283,2.78973693455129,2.62745564169631,46.10583195234,43.1607234896162,41.1417639432853,40.6162182241421,39.7789855858654,36.1525839369613,32.8861031569682,31.7655484997385,27.6481564270394,25.1197654350362,24.9807758798591,23.5622458758418,22.3917397864646,21.900021891676,21.7646852982626,21.6704789610659,20.8935912726363,20.78148179504,20.3081649787342,20.1367875165235,19.8478949263734,19.439099730812,18.9362702957842,18.3579029537197,18.1706562305444,15.854544826175,15.67883572726,14.6180137170381,14.3311062863934,13.5344717127731,12.9215645955221,11.8688491564598,11.1578645158304,10.8453435829731,10.1766043347684,9.34806928051863,8.522994826618,8.49288755532288,8.05855543024414,7.96008384905675,7.73567684844406,5.25534712702971,4.8113784889159,4.6397084513437,4.25750644345912,3.50191990208325,2.92935888911263,2.68743040370355,21.1729697316108,9.42328666351578,3.70079401432848,3.44540646806805,2.35225867515329,2.25976708063955,5.70521680966971,3.20427492902515,60.0362308806446,18.6015942045366,17.8639342436323,17.8434073714439,17.8126760727161,17.7950558141117,17.4545924703922,17.3071465749948,17.1936703619909,17.1681847748447,17.1186819642422,17.1116363453666,17.1102274030416,17.0332503959747,16.9763730444757,16.9537900972103,16.9171437492488,16.9019508067597,16.8623566324271,16.8217635160776,16.8175045483819,16.8050799566303,16.7695079288439,16.6942391737494,16.6564973735799,16.6502461322852,16.5656278687299,16.513720103301,16.5080500099678,16.4899401811273,16.3987198181955,16.3508007043865,16.2713040538615,16.1493281434008,16.1422202730942,16.0660502124752,16.0573258902284,15.987892400769,15.9352416438171,15.9176663288293,15.7064076576521,15.5468375298472,15.5466365062712,15.5359842016412,15.4400849516219,15.4339430193406,14.9876530007503,14.9636693680233,22.8567432425572,22.8040258781381,22.7752047640768,22.600753100535,21.26084718376,21.0477373761922,20.941383288413,20.2692149084912,19.9698646934514,19.8739349560792,19.6446304322614,19.4820400187451,18.6986827089862,17.6281539071434,16.0812073773676,15.9435573712599,13.9857104931239,9.31161927723081,4.06011169260699,3.08725473489223,4.82851412413945,3.62338474977333,2.7182923444479,4.82851412413945,3.62338474977333,2.7182923444479,3.89300819230754,3.38170000646001,2.87381372541297,10.6370456615209,9.16947281168527,8.1618431988927,8.09827845050814,7.96505221424297,5.717099398135,3.32180630047215,16.7769867401313,3.7638280318425,3.51748478321809,3.46056427117438,2.50663986351731,2.33359184059253,22.9713556444982,22.7030741927304,8.96427736026615,8.70552051597353,12.5924794344816,11.9814215696202,10.8698804675621,10.1771522960867,7.93200317859503,7.81331144540585,4.71428189158595,4.14078220474963,3.69854391885719,3.48781329054908,2.44502954491056,2.29601777557942,2.240925083768,12.5924794344816,11.9814215696202,10.8698804675621,10.1771522960867,7.93200317859503,7.81331144540585,4.71428189158595,4.14078220474963,3.69854391885719,3.48781329054908,2.44502954491056,2.29601777557942,2.240925083768,30.8732891020125,6.59228748517027,4.52796094012902,2.68735723600618,2.29274023213369,2.35391368269323,2.35391368269323,2.22204068529081,2.21216723311021,2.16217287414067]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> \u003c/th>\n      <th>GSE\u003c/th>\n      <th>GSM\u003c/th>\n      <th>GPL\u003c/th>\n      <th>TGROUP\u003c/th>\n      <th>FOUND_TGROUP\u003c/th>\n      <th>ORGANISM\u003c/th>\n      <th>EXP_SCORE\u003c/th>\n      <th>FOUND_SCORE\u003c/th>\n      <th>score_ratio\u003c/th>\n    \u003c/tr>\n  \u003c/thead>\n\u003c/table>","options":{"columnDefs":[{"className":"dt-right","targets":[7,8,9]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->



### Using the bioqc_solid signature set
We can re-run the above analysis with the signatures provided by the BioQC authors. These signatures are not as toroughly validated as the 'gtex solid' set used above. Nonetheless it is interesting (1) to compare the results and (2) to take additional tissues into account. 


```r
sql = "
select /*+ parallel(16) */ *
from bioqc_contamination
where tissue_set = 'bioqc_solid'
"
res_bqc = data.table(dbGetQuery(mydb, sql))
res_bqc = res_bqc[,MIN_FOUND_PVALUE:=as.numeric(MIN_FOUND_PVALUE)]
res_bqc = res_bqc[,MIN_EXP_PVALUE:=as.numeric(MIN_EXP_PVALUE)]
res_bqc = res_bqc[,EXP_SCORE:=absLog10p(MIN_EXP_PVALUE)]
res_bqc = res_bqc[,FOUND_SCORE:=absLog10p(MIN_FOUND_PVALUE)]
res_bqc = res_bqc[,total_tissue:=length(unique(GSM)), by=c("TGROUP")]
#res[,cnt_per_group:=length(unique(GSM)), by=c("GPL", "TGROUP")]
#res = res[cnt_per_group >= 100]
pbonf_bqc = 0.05 / nrow(res_bqc)
```

Again, we note that the signatures in general are highly specific.

```r
ggplot(res_bqc, aes(x=FOUND_TGROUP, y=FOUND_SCORE)) + 
  geom_boxplot() + facet_wrap(~TGROUP) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-10-1.png" width="1152" style="display:block; margin: auto" style="display: block; margin: auto;" />

We filter the results according to the cutoff

```r
res_bqc = res_bqc[,score_ratio:=FOUND_SCORE - EXP_SCORE]
fil_bqc = res_bqc[TGROUP != FOUND_TGROUP & score_ratio > 2]

#### exclude 'double contaminations'
fil_bqc = fil_bqc[, rk:=frankv(FOUND_SCORE, order=-1), by=c("GSM", "GPL", "TGROUP")]
fil_bqc = fil_bqc[rk==1]
```

and display the 'contamination matrix' 

```r
aggr = sqldf('select TGROUP, FOUND_TGROUP, count(GSM) as cnt from fil_bqc group by TGROUP, FOUND_TGROUP')
ggplot(aggr, aes(y=TGROUP, x=FOUND_TGROUP)) +
  geom_tile(aes(fill=cnt)) +
  geom_text(aes(label=cnt))
```

<img src="50_contamination_results_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display:block; margin: auto" style="display: block; margin: auto;" />


Investigtation why this approach shows many more hits for blood and muscle contamination 
has shown that the positive hits are mainly due to single-cell blood signatures (e.g. Monocytes) and the 'smooth muscle' muscle signature which has not been profiled in the gtex dataset (data not shown).



## signature-free method?
To find more subtle contamination, we could conduct an exemplary study searching for pancreas contamination in kidney. This is a likely contamination in mouse due to the spatial spatial closeness of the organs, as already noted by the BioQC authors. We could do so by assessing the abundance of the Elastase and Amylase genes which are only expressed in pancreas.
