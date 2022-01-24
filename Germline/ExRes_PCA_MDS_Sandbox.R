---
title: "ExRes PCA/MDS sandbox"
output: html_notebook
---

### 7. nMDS-ExAC Z scores: A bundled measure of functional constraints per gene

The ExAC Z scores (read more at where?) are:

> ExAC_delZ,	ExAC_dupZ,	ExAC_cnvZ,	ExAC_synZ,	ExAC_misZ


Since they are Z scores ... they should approximately follow known distribution (SND)...right?


#### Workflow

--> Extract all ExAC data into numeric mtx
```{r}
# selecting a subset of SV size < 20 kb 
# as we are evaluating whether ExAC Z scores together are a good representation/estimate of 
# "functional constraint" corresponding to a gene-variant pair. We are hence omitting larger variants
# that likely span multiple genes.

# Also, vegan::metaMDS doesn't seem tolerate huge data frames well.

elite_variants_exres %>%
  subset(select=c("ExAC_delZ",	"ExAC_dupZ",	"ExAC_cnvZ",	"ExAC_synZ",	"ExAC_misZ")) -> Z.mtx
```

--> Already standardised

--> Filter out NA, beware of zeros
```{r}
Z.mtx[complete.cases(Z.mtx), ] -> Z.mtx ## removes all rows containing any NA

round(Z.mtx$ExAC_cnvZ, 9) -> Z.mtx$ExAC_cnvZ
```

```{r}
# Subsample to avoid overtime issue
Z.mtx[sample(nrow(Z.mtx), 400), ] -> Z.mtx.sample
```

--> Run MDS compression
```{r}
MDS.S <- vegan::metaMDS(Z.mtx.sample, distance = "euclidean", autotransform = FALSE, trace = FALSE)
```


--> Check stress level of MDS
```{r}
MDS.S$stress ## Checking MDS stress for interpretation
MDS.S$stress <= 0.20
```


Recall the standard level of acceptance of MDS stress = 0.20.


--> Obtain Cartesian coordinate in S_mtx
```{r}
MDS.xy <- data.frame(MDS.S$points) ## Sending the cartesian coordinates in MDS profile to MDS.xy
```


--> Send coordinate to new dummy df + assign factors of interest (e.g. Caller)

#### Assigning factors: Great code
```{r}
# precious codes
MDS.xy$Caller = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"Caller"]
names(MDS.xy)[names(MDS.xy) == 'Caller.Caller'] <- 'Caller'

MDS.xy$GnomAD_pLI = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"GnomAD_pLI"]
names(MDS.xy)[names(MDS.xy) == 'GnomAD_pLI.GnomAD_pLI'] <- 'GnomAD_pLI'
MDS.xy$ExAC_pLI = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"ExAC_pLI"]
names(MDS.xy)[names(MDS.xy) == 'ExAC_pLI.ExAC_pLI'] <- 'ExAC_pLI'

MDS.xy$SV_type = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"SV_type"]
names(MDS.xy)[names(MDS.xy) == 'SV_type.SV_type'] <- 'SV_type'

MDS.xy$SV_length = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"SV_length"]
names(MDS.xy)[names(MDS.xy) == 'SV_length.SV_length'] <- 'SV_length'

MDS.xy$SV_chrom = elite_variants_exres[match(row.names(MDS.xy), row.names(elite_variants_exres)),"SV_chrom"]
names(MDS.xy)[names(MDS.xy) == 'SV_chrom.SV_chrom'] <- 'SV_chrom'
```


--> Run MDS plot and annotate with factor levels 

#### 7-1. ExACZ ~ Caller
```{r}
# MDS: (1) By Caller
ggplot(MDS.xy, aes(MDS1, MDS2, color = Caller)) + geom_point(size = 1) +

  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "grey15", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey30", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 13.5),
        legend.key = element_rect(fill = NA)) +
  
  guides(color = guide_legend(override.aes = list(size = 4.5)))
```



#### 7-2. ExACZ ~ ACMG class
```{r}
# MDS: (2) By ACMG class of SV
MDS.xy$ACMG_class <- as.factor(MDS.xy$ACMG_class) ## coerce as factor level (discrete) instead of numeric (continuous)

ggplot(MDS.xy, aes(MDS1, MDS2, color = ACMG_class)) + geom_point() + theme_bw() +
  
  scale_colour_viridis(discrete = TRUE, name = "ACMG class") +
  
    theme(legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "grey97", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey86", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) + 
  
  guides(color = guide_legend(override.aes = list(size = 4.5)))
```


#### 7-3. ExACZ ~ SV type
```{r}
# MDS: (3) By SV type
ggplot(MDS.xy, aes(MDS1, MDS2, color = SV_type)) + geom_point(size = 1) +
  
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "grey12", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey30", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA)) +
  
  scale_color_discrete(name = "SV type") +
  
  guides(color = guide_legend(override.aes = list(size = 4.5)))
```


```{r}
# Try facet_wrap
MDS.xy %>% filter(SV_type != "SVA") %>%
  ggplot(.) + aes(MDS1, MDS2, color = SV_type) + geom_point(size = 1) +
  
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "grey12", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey30", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA)) +
  
  guides(color = guide_legend(override.aes = list(size = 4.5))) +
  
  scale_color_discrete(name = "SV type") +
  
  facet_wrap( ~SV_type)
```



#### 7-4. ExACZ ~ SV length
```{r}
# MDS: (4) By SV length
MDS.xy %>% filter(!is.na(SV_length)) %>% mutate(SV_length = log10(abs(SV_length))) %>%

ggplot(aes(MDS1, MDS2, color = SV_length)) + geom_point(size = 1) +
  
  scale_color_viridis(option = "A") +
  
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "grey12", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey30", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA)) +
  
  guides(color = guide_legend(override.aes = list(size = 4.5)))
```


#### 7-5. ExACZ ~ variant location
```{r}
# MDS: (5) By variant location
ggplot(MDS.xy, aes(MDS1, MDS2, color=SV_chrom)) + geom_point(size = 1) +
  
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "grey12", colour = "black", size = 1),
        #panel.grid.major = element_line(color = "grey30", size=0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA)) +
  
  scale_color_discrete(name = "Chromosome") +
  
  guides(color = guide_legend(override.aes = list(size = 4.5)))
```


### 8. PCA:

#### 8-1. ExAC_Z continued: ACMG score and ExAC Z measure of functional constraint
```{r}
# using entire data frame instead of subsample as in MDS restrictions
Z.pca <- princomp(Z.mtx)
```

```{r}
# first-stage analysis
plot(Z.pca$scores) ## shows visuals of what first two PCs look like. can specify which = c(a,b)
screeplot(Z.pca) ## %var explained by each PC
summary(Z.pca) ## important statistics of each PC, mostly %var accounted
loadings(Z.pca) ## "loadings" in PCA refers to the linear (de)composition of each PC
```

Select PC1,2 = ~ 80%; PC1,2,3 = ~ 91% (Quite low!)

#### 2D PCA (80%)
```{r}
# PCA plot: 2D
Z.coor <- data.frame(Z.pca$scores)

# assigning factors
Z.coor$Caller = large.full[match(row.names(Z.coor), row.names(large.full)),"Caller"]
Z.coor$ACMG_class = large.full[match(row.names(Z.coor), row.names(large.full)),"ACMG_class"]
Z.coor$SV_type = large.full[match(row.names(Z.coor), row.names(large.full)),"SV_type"]
```


```{r}
# Plot 2D PCA
## by SV type
ggplot(Z.coor) + aes(Comp.1, Comp.2, color=SV_type) + geom_point() +
   theme_bw()

## by Caller
ggplot(Z.coor) + aes(Comp.1, Comp.2, color=Caller) + geom_point() +
   theme_bw()

## by ACMG class
Z.coor$ACMG_class <- as.factor(Z.coor$ACMG_class) ## discrete leveling

ggplot(Z.coor) + aes(Comp.1, Comp.2, color=ACMG_class) + geom_point() +
  scale_colour_viridis(discrete = TRUE) +
   theme_bw()
```


Consider subsample to n=600:

#### Subsampling: n=600
```{r}
# sample down to n=600
Z.coor[sample(nrow(Z.coor), 600), ] -> Z.pca.sample

# Plot 2D PCA
## by SV type
ggplot(Z.pca.sample) + aes(Comp.1, Comp.2, color=SV_type) + geom_point() +
   theme_bw()

## by Caller
ggplot(Z.pca.sample) + aes(Comp.1, Comp.2, color=Caller) + geom_point() +
   theme_bw()

## by ACMG class
ggplot(Z.pca.sample) + aes(Comp.1, Comp.2, color=ACMG_class) + geom_point() +
  scale_colour_viridis(discrete = TRUE) +
   theme_bw()
```


```{r}
Z.pca.sample %>% filter(SV_type != "SVA") %>%
  ggplot(.) + aes(Comp.1, Comp.2, color=SV_type) + geom_point() +
   theme_bw() + facet_wrap( ~SV_type)
```






    











