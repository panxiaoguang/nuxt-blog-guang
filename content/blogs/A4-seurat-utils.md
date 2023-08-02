---
title: Seurat3.1的灵活操作指南
date: 20th Nov 2019
description: Seurat3.1的灵活操作指南.
image: /blogs-img/singlecell.webp
alt: Seurat3.1的灵活操作指南
ogImage: /blogs-img/singlecell.webp
tags: ['单细胞测序', '生物信息', 'r语言']
published: true
---

官网3.1版本已经无法找到该指南的链接，其实还是有的，网址：
[https://satijalab.org/seurat/v3.1/interaction_vignette.html](https://satijalab.org/seurat/v3.1/interaction_vignette.html)

### 载入数据

下面演示了一些与Seurat对象进行交互的有用功能。出于演示目的，我们将使用在第一个指导教程中创建的2700 PBMC对象。您可以在此处下载预先计算的对象。为了模拟有两个重复的情况，将一半命名为“rep1",另一半命名为"rep2"

```r
library(Seurat)
pbmc <- readRDS(file = "../data/pbmc3k_final.rds")

### 随机设置两个重复
set.seed(42)
pbmc$replicate <- sample(c("rep1", "rep2"), size = ncol(pbmc), replace = TRUE)
```
### 从细胞聚类和样本重复中切换Idents

```r
# 默认画的是object@ident)
DimPlot(pbmc, reduction = "umap")
```
![](https://cdn.jsdelivr.net/gh/panxiaoguang/MyImage/img/exam_dim.png)
```r
# 把细胞分类因子先储存到对象中
pbmc$CellType <- Idents(pbmc)
# 切换Idents
Idents(pbmc) <- "replicate"
DimPlot(pbmc, reduction = "umap")
```
![](https://cdn.jsdelivr.net/gh/panxiaoguang/MyImage/img/sample_class.png)

```r
# alternately : DimPlot(pbmc, reduction = 'umap', group.by = 'replicate') you can pass the
# shape.by to label points by both replicate and cell type

# Switch back to cell type labels
Idents(pbmc) <- "CellType"
```

### 分别或同时统计不同聚类或者不同样本来源的细胞数目
```r
# 每个聚类包含多少细胞？
table(Idents(pbmc))
```

```
## 
##  Naive CD4 T Memory CD4 T   CD14+ Mono            B        CD8 T 
##          697          483          480          344          271 
## FCGR3A+ Mono           NK           DC     Platelet 
##          162          155           32           14
```
```r
# 每组重复包含多少细胞？
table(pbmc$replicate)
```

```
## 
## rep1 rep2 
## 1348 1290
```
```r
# 每个聚类细胞数占比
prop.table(table(Idents(pbmc)))
```
```
## 
##  Naive CD4 T Memory CD4 T   CD14+ Mono            B        CD8 T 
##  0.264215315  0.183093252  0.181956027  0.130401820  0.102729340 
## FCGR3A+ Mono           NK           DC     Platelet 
##  0.061410159  0.058756634  0.012130402  0.005307051
```
```r
# 样本分组和细胞聚类一起统计
table(Idents(pbmc), pbmc$replicate)
```
```
##               
##                rep1 rep2
##   Naive CD4 T   354  343
##   Memory CD4 T  249  234
##   CD14+ Mono    232  248
##   B             173  171
##   CD8 T         154  117
##   FCGR3A+ Mono   81   81
##   NK             81   74
##   DC             18   14
##   Platelet        6    8
```
```r
prop.table(table(Idents(pbmc), pbmc$replicate), margin = 2)
```
```
##               
##                       rep1        rep2
##   Naive CD4 T  0.262611276 0.265891473
##   Memory CD4 T 0.184718101 0.181395349
##   CD14+ Mono   0.172106825 0.192248062
##   B            0.128338279 0.132558140
##   CD8 T        0.114243323 0.090697674
##   FCGR3A+ Mono 0.060089021 0.062790698
##   NK           0.060089021 0.057364341
##   DC           0.013353116 0.010852713
##   Platelet     0.004451039 0.006201550
```

### 提取特定的Seurat子集做亚型分析
```r
# What are the cell names of all NK cells?
WhichCells(pbmc, idents = "NK")
```
```
##   [1] "AAACCGTGTATGCG" "AAATTCGATTCTCA" "AACCTTACGCGAGA" "AACGCCCTCGTACA"
##   [5] "AACGTCGAGTATCG" "AAGATTACCTCAAG" "AAGCAAGAGCTTAG" "AAGCAAGAGGTGTT"
##   [9] "AAGTAGGATACAGC" "AATACTGAATTGGC" "AATCCTTGGTGAGG" "AATCTCTGCTTTAC"
##  [13] "ACAAATTGTTGCGA" "ACAACCGAGGGATG" "ACAATTGATGACTG" "ACACCCTGGTGTTG"
##  [17] "ACAGGTACTGGTGT" "ACCTGGCTAAGTAG" "ACGAACACCTTGTT" "ACGATCGAGGACTT"
##  [21] "ACGCAATGGTTCAG" "ACGCTGCTGTTCTT" "ACGGAACTCAGATC" "ACGTGATGTGACAC"
##  [25] "ACGTTGGAGCCAAT" "ACTGCCACTCCGTC" "ACTGGCCTTCAGTG" "ACTTCAACGTAGGG"
##  [29] "AGAACAGAAATGCC" "AGATATACCCGTAA" "AGATTCCTGTTCAG" "AGCCTCTGCCAATG"
##  [33] "AGCGATTGAGATCC" "AGGATGCTTTAGGC" "AGGGACGAGTCAAC" "AGTAATACATCACG"
##  [37] "AGTCACGATGAGCT" "AGTTTGCTACTGGT" "ATACCACTGCCAAT" "ATACTCTGGTATGC"
##  [41] "ATCCCGTGCAGTCA" "ATCTTTCTTGTCCC" "ATGAAGGACTTGCC" "ATGATAACTTCACT"
##  [45] "ATGATATGGTGCTA" "ATGGACACGCATCA" "ATGGGTACATCGGT" "ATTAACGATGAGAA"
##  [49] "ATTCCAACTTAGGC" "CAAGGTTGTCTGGA" "CAATCTACTGACTG" "CACCACTGGCGAAG"
##  [53] "CACGGGTGGAGGAC" "CAGATGACATTCTC" "CAGCAATGGAGGGT" "CAGCGGACCTTTAC"
##  [57] "CAGCTCTGTGTGGT" "CAGTTTACACACGT" "CATCAGGACTTCCG" "CATCAGGATAGCCA"
##  [61] "CATGAGACGTTGAC" "CATTACACCAACTG" "CATTTCGAGATACC" "CCTCGAACACTTTC"
##  [65] "CGACCACTAAAGTG" "CGACCACTGCCAAT" "CGAGGCTGACGCTA" "CGCCGAGAGCTTAG"
##  [69] "CGGCGAACGACAAA" "CGGCGAACTACTTC" "CGGGCATGTCTCTA" "CGTACCTGGCATCA"
##  [73] "CGTGTAGACGATAC" "CGTGTAGAGTTACG" "CGTGTAGATTCGGA" "CTAAACCTCTGACA"
##  [77] "CTAACGGAACCGAT" "CTACGCACTGGTCA" "CTACTCCTATGTCG" "CTAGTTACGAAACA"
##  [81] "CTATACTGCTACGA" "CTATACTGTCTCAT" "CTCGACTGGTTGAC" "CTGAGAACGTAAAG"
##  [85] "CTTTAGTGACGGGA" "GAACCAACTTCCGC" "GAAGTGCTAAACGA" "GAATGCACCTTCGC"
##  [89] "GAATTAACGTCGTA" "GACGGCACACGGGA" "GAGCGCTGAAGATG" "GAGGTACTGACACT"
##  [93] "GAGGTGGATCCTCG" "GATAGAGAAGGGTG" "GATCCCTGACCTTT" "GCACACCTGTGCTA"
##  [97] "GCACCACTTCCTTA" "GCACTAGAGTCGTA" "GCAGGGCTATCGAC" "GCCGGAACGTTCTT"
## [101] "GCCTACACAGTTCG" "GCGCATCTTGCTCC" "GCGCGATGGTGCAT" "GGAAGGTGGCGAGA"
## [105] "GGACGCTGTCCTCG" "GGAGGCCTCGTTGA" "GGCAAGGAAAAAGC" "GGCATATGCTTATC"
## [109] "GGCCGAACTCTAGG" "GGCTAAACACCTGA" "GGGTTAACGTGCAT" "GGTGGAGAAACGGG"
## [113] "GTAGTGTGAGCGGA" "GTCGACCTGAATGA" "GTGATTCTGGTTCA" "GTGTATCTAGTAGA"
## [117] "GTTAAAACCGAGAG" "GTTCAACTGGGACA" "GTTGACGATATCGG" "TAACTCACTCTACT"
## [121] "TAAGAGGACTTGTT" "TAATGCCTCGTCTC" "TACGGCCTGGGACA" "TACTACTGATGTCG"
## [125] "TACTCTGAATCGAC" "TACTGTTGAGGCGA" "TAGCATCTCAGCTA" "TAGCCCACAGCTAC"
## [129] "TAGGGACTGAACTC" "TAGTGGTGAAGTGA" "TAGTTAGAACCACA" "TATGAATGGAGGAC"
## [133] "TATGGGTGCATCAG" "TATTTCCTGGAGGT" "TCAACACTGTTTGG" "TCAGACGACGTTAG"
## [137] "TCCCGAACACAGTC" "TCCTAAACCGCATA" "TCGATTTGCAGCTA" "TCTAACACCAGTTG"
## [141] "TGATAAACTCCGTC" "TGCACAGACGACAT" "TGCCACTGCGATAC" "TGCTGAGAGAGCAG"
## [145] "TGGAACACAAACAG" "TGGTAGACCCTCAC" "TGTAATGACACAAC" "TGTAATGAGGTAAA"
## [149] "TTACTCGATCTACT" "TTAGTCTGCCAACA" "TTCCAAACTCCCAC" "TTCCCACTTGAGGG"
## [153] "TTCTAGTGGAGAGC" "TTCTGATGGAGACG" "TTGTCATGGACGGA"
```
```r
# 提取NK细胞的表达矩阵
nk.raw.data <- as.matrix(GetAssayData(pbmc, slot = "counts")[, WhichCells(pbmc, ident = "NK")])

# 获取基因表达量大于1的对象
subset(pbmc, subset = MS4A1 > 1)
```
```
## An object of class Seurat 
## 13714 features across 414 samples within 1 assay 
## Active assay: RNA (13714 features)
##  2 dimensional reductions calculated: pca, umap
```
```r
subset(pbmc, subset = replicate == "rep2")
```
```
## An object of class Seurat 
## 13714 features across 1290 samples within 1 assay 
## Active assay: RNA (13714 features)
##  2 dimensional reductions calculated: pca, umap
```
```r
# 选择两个细胞类型
subset(pbmc, idents = c("NK", "B"))
```
```r
## An object of class Seurat 
## 13714 features across 499 samples within 1 assay 
## Active assay: RNA (13714 features)
##  2 dimensional reductions calculated: pca, umap
```
```r
# 排除掉某些细胞类型
subset(pbmc, idents = c("NK", "B"), invert = TRUE)
```
```
## An object of class Seurat 
## 13714 features across 2139 samples within 1 assay 
## Active assay: RNA (13714 features)
##  2 dimensional reductions calculated: pca, umap
```
```
# note that if you wish to perform additional rounds of clustering after subsetting we recommend
# re-running FindVariableFeatures() and ScaleData()
```

### 计算基因的平均表达量
```r
# 计算平均表达量
cluster.averages <- AverageExpression(pbmc)
head(cluster.averages[["RNA"]][, 1:5])
```

|  | Native CD4 T | Memory CD4 T | CD14+Mono | B | CD8 T |
| --- | --- | --- | --- | --- | --- |
| AL627309.1<span class="Apple-tab-span" style="white-space:pre"></span> | 0.0061287 | 0.0059273 | 0.0485434 | 0.0000000 | 0.0205459 |
| AP006222.2 | 0.0000000 | 0.0082061 | 0.0108847 | 0.0000000 | 0.0119149 |
| RP11-206L10.2 | 0.0074531 | 0.0000000 | 0.0000000 | 0.0206503 | 0.0000000 |
| RP11-206L10.9 | 0.0000000 | 0.0000000 | 0.0105012 | 0.0000000 | 0.0000000 |
| LINC00115 | 0.0191189 | 0.0246905 | 0.0375374 | 0.0388854 | 0.0194828 |
| NOC2L | 0.4974632 | 0.3598115 | 0.2725375 | 0.5865349 | 0.5570490 |


```r
# 返回Seurat对象用于下游分析
orig.levels <- levels(pbmc)
Idents(pbmc) <- gsub(pattern = " ", replacement = "_", x = Idents(pbmc))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(pbmc) <- orig.levels
cluster.averages <- AverageExpression(pbmc, return.seurat = TRUE)
cluster.averages
```
```
## An object of class Seurat 
## 13714 features across 9 samples within 1 assay 
## Active assay: RNA (13714 features)
```
```r
# How can I plot the average expression of NK cells vs. CD8 T cells?  Pass do.hover = T for an
# interactive plot to identify gene outliers
CellScatter(cluster.averages, cell1 = "NK", cell2 = "CD8_T")
```
![](https://cdn.jsdelivr.net/gh/panxiaoguang/MyImage/img/20191206152652.png)

```r
# How can I calculate expression averages separately for each replicate?
cluster.averages <- AverageExpression(pbmc, return.seurat = TRUE, add.ident = "replicate")
CellScatter(cluster.averages, cell1 = "CD8_T_rep1", cell2 = "CD8_T_rep2")
```
![](https://cdn.jsdelivr.net/gh/panxiaoguang/MyImage/img/20191206152710.png)

```r
# You can also plot heatmaps of these 'in silico' bulk datasets to visualize agreement between
# replicates
DoHeatmap(cluster.averages, features = unlist(TopFeatures(pbmc[["pca"]], balanced = TRUE)), size = 3, 
    draw.lines = FALSE)
```

![](https://cdn.jsdelivr.net/gh/panxiaoguang/MyImage/img/20191206152720.png)