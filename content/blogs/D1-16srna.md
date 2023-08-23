---
title: 16S rRNA基因扩增子分析 
date: 23th Aug 2023
description: 在此工作流程中，介绍了 Qiime2 和 R 中 16S rRNA 基因扩增子数据分析的主要步骤。本教程是为哥本哈根大学食品科学系的 MAC 2023 课程准备的。尽管这些步骤是为 Oxford Nanopore Tech (ONT) 测序设计的，但也在 Ilumina 短读长上进行了测试。
image: https://picshack.net/ib/qlbONcTVD5.png
alt: 16s rRNA sequencing
ogImage: https://picshack.net/ib/qlbONcTVD5.png
tags: ['生物信息','R语言']
published: true
---

在此工作流程中，介绍了 Qiime2 和 R 中 16S rRNA 基因扩增子数据分析的主要步骤。本教程是为哥本哈根大学食品科学系的 MAC 2023 课程准备的。尽管这些步骤是为 Oxford Nanopore Tech (ONT) 测序设计的，但也在 Ilumina 短读长上进行了测试。

### Qiime2中的步骤

#### 1. 制作清单文件。

该文件对于将原始读数导入 Qiime2 非常重要。该文件只是一个制表符分隔的文本，有两列（如果我们有双端读取）；正向绝对文件路径和反向绝对文件路径。每列包含正向和反向通道每次读取的绝对路径。

```bash
#path to forward
ls -f ~/microbiome_analysis_ku/data/*_R1_001.fastq.gz | sort -V > r1path
#path to reverse
ls -f ~/microbiome_analysis_ku/data/*_R2_001.fastq.gz | sort -V > r2path

#checking if we have all pairs of forward and reverse reads
diff -s <(cat r1path | cut -d_ -f4) <(cat r2path | cut -d_ -f4) #if the outcome is "identical", then we good :)

#Cutting sample IDs from the path file (from either is good)
cut -d_ -f4 r1path > SampleID

#adding sample ids to paths
paste SampleID r1path r2path > manifest

#Binning ids and forward reverse columns together with tab-delimited seperation
sed -i $'1 i\\\nsampleid \t forward-absolute-filepath \t reverse-absolute-filepath' manifest
```

一些数据的推理调查

```bash
## Exploring the reads
zcat ./data/MAC2023_iSeq001_S1_L001_R2_001.fastq.gz 


zcat ./data/MAC2023_iSeq001_S1_L001_R2_001.fastq.gz | grep ^@FS10000714| cut -d ":" -f1 | wc -l

wc -l ./data/*_R1_001.fastq.gz | awk '{$1=$1};1'| cut -d" " -f1|datamash min 1
wc -l ./data/*_R1_001.fastq.gz | awk '{$1=$1};1'| cut -d" " -f1|head -n -1|datamash max 1
wc -l ./data/*_R1_001.fastq.gz | awk '{$1=$1};1'| cut -d" " -f1|head -n -1|datamash mean 1

#reverse
wc -l ./data/*_R2_001.fastq.gz| awk '{$1=$1};1' | cut -d " " -f1|datamash min 1
wc -l ./data/*_R2_001.fastq.gz| awk '{$1=$1};1' | cut -d " " -f1|datamash max 1
wc -l ./data/*_R2_001.fastq.gz| awk '{$1=$1};1' | cut -d " " -f1|datamash mean 1

#searching for primers

zcat ./data/MAC2023_iSeq001_S1_L001_R2_001.fastq.gz | head
zgrep '^CCTACGGG.GGC.GCAG' ./data/MAC2022_iSeq001_S1_L001_R2_001.fastq.gz | wc -l
zgrep '^GACTAC..GGGTATCTAATCC' ./data/MAC2022_iSeq001_S1_L001_R2_001.fastq.gz | wc -l

#Eventually you can use Seqkit package

seqkit stats ./data/MAC2023_iSeq001_S1_L001_R2_001.fastq.gz
```

#### 2. 使用创建的清单文件将读取导入 Qiime2

```bash
source activate qiime2.X
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path ~/MAC2023/manifest.tsv \                    # link to the folder in which my manifest file is located
  --output-path ~/MAC2023/demuxed_MAC23.qza                 # link to the path where I want my demultiplexed data to be exported in

```

#### 3. 可视化解复用读取的推论统计

```bash
qiime demux summarize \
  --i-data ./demuxed_mac23.qza \
  --o-visualization ./demuxed_mac23.qzv
```

#### 4. qiime2中DADA2包的过滤、去重复、样本推断、嵌合体识别、双端读段合并。

```bash
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ~/demuxed_mac23.qza \         #the input file for denoising is our demultiplexed pairedEnd reads that was generated in the previous step.
--p-trim-left-f 17 \                                        #length of my forward primer (17 nt) 
--p-trim-left-r 21 \                                        #length of my reverse primer (21 nt)
--p-trunc-len-f 260 \                                       #truncation length for forward reads. I.e. we truncate reads over 260 base since their quality started dropping from this point onwards. 
--p-trunc-len-r 220 \                                       #truncation length for reverse reads. I.e. we truncate reads over 220 base since their quality started dropping from this point onwards. 
--o-table ~/count_table.qza \                      #the ASV count table
--o-representative-sequences ~/repseqs.qza \ #the representative sequences for in each read
--o-denoising-stats ~/denoising_stats.qza \  #the status of the denoising process in a full catalogue
--p-n-threads 10
```

##### 4.1. 可视化计数表

```bash
qiime feature-table summarize --i-table ~/count_table.qza \ #The ASV table as the input
--m-sample-metadata-file ~/metadataA.tsv \           #The metadata as the input
--o-visualization ./table.qzv                         #The qzv format of our ASV table
```

#### 4.2. 可视化代表性序列

```bash
qiime feature-table tabulate-seqs \
--i-data ~/repseqs.qza \
--o-visualization ~/repseqs.qzv
```

#### 4.3. 可视化去噪统计数据

```bash
qiime metadata tabulate \
--m-input-file ./denoising_stats.qza \
--o-visualization ./denoising_stats.qzv
```

#### 5. 通过朴素贝叶斯方法训练基于引物的区域特异性分类器（在 Qiime2 中）

```bash
qiime feature-classifier fit-classifier-naive-bayes \   # here you can train your classifier based on a Naive-Bayes method
--i-reference-reads ~/derepseqs-uniq-341f-805r.qza \    # Dereplicated sequences based on the primer set as input
--i-reference-taxonomy ~/dereptaxa-uniq-341f-805r.qza\  # Dereplicated taxonomic annotations based on the primer set as input
 --o-classifier ~/silva-classifier-primered4.qza  

# Running the classifier
qiime feature-classifier classify-sklearn \                               # Sklearn package for classification
--i-reads ~/repseqs.qza \                                        # Representative sequences as the (input)
--i-classifier ~/classifier/silva138-classifier-341f-805r.qza \  # Our costumized SILVA 138 classifier (input)
--o-classification ~/taxonomy.qza \                     # The taxonomic annotation linked to the repseqs (output)
--p-n-jobs 10                                                             # The number of cores/jobs
```

##### 5.1. 可视化分类表

```bash
qiime metadata tabulate \
--m-input-file ~/Taxonomy/taxonomy.qza \
--o-visualization ~/Taxonomy/taxonomy.qzv
```

#### 6. 使用支持 SATE 的系统发育放置 (SEPP) 方法创建系统发育树以进行多样性分析

```bash
qiime fragment-insertion sepp \                             # The package for generating the phylogenetic tree
--i-representative-sequences ~/repseqs.qza \                # Our representative sequences (input)
--i-reference-database ~/epp-ref-gg-13-8.qza \              # The SEPP sequence dataset (input)
--o-tree ~/tree.qza \                                       # The rooted tree (output)
--o-placements ~/dss/tree-placements.qza \
--p-threads 10  
```

### R 中的步骤

#### 1.安装和加载库

```r
#Update R
install.packages('installr')
installr::updateR()

github.pkg <- c("jfukuyama/phyloseqGraphTest", "jbisanz/qiime2R") 
bioc.pkg <- c("phyloseq", "DESeq2", "MicrobiotaProcess", "ggtree", "vsn")
cran.pkg <- c("tidyverse", "readr", "ape", "pacman", "picante", "glue", "vegan", "devtools", "ggrepel", "reshape2", "BiocManager",
              "ggnetwork", "DT", "VennDiagram", "lsmeans", "pheatmap", "phyloseqGraphTest")

inst.pkg <- cran.pkg %in% installed.packages()

BiocManager::install("ggtree")

if (any(!inst.pkg)){ install.packages(cran.pkg[!inst.pkg],repos = "http://cran.rstudio.com/") } 
 
inst.pkg <- github.pkg %in% installed.packages() 
if (any(!inst.pkg)){ devtools::install_github(github.pkg[!inst.pkg], force = TRUE) } 


inst.pkg <- bioc.pkg %in% installed.packages() 
if(any(!inst.pkg)){ BiocManager::install(bioc.pkg[!inst.pkg])
}
```

```r
pkg = c("tidyverse", "phyloseq", "DESeq2", "qiime2R" , "gridExtra", "BiocManager", 
"vegan", "ggtree", "devtools", "ggrepel", "reshape2", "ggnetwork", "igraph", "biomformat", 
"pheatmap", "glue", "ape", "readr", "vsn", "vegan"
)

for(i in 1:length(pkg)){

  library(pkg[i], character.only = TRUE, verbose = FALSE, attach.required = FALSE)

}
```

```r
save.image('./data/r.RData')
load("./data/r.RData")
```

#### 2.从qiime2导入工件到r
```r
#making phyloseq objects from qiime files
ps <- qiime2R::qza_to_phyloseq(features = "~/data/ccd/table-ccd.qza", 
                       taxonomy = "~/data/ccd/taxonomy-ccd.qza", 
                     tree = "./tree-ccd.qza")
repseqs <- qiime2R::read_qza("~/data/ccd/repseqa.qza")$data
#we need to merge the refseqs like this since in above fun, there is no argument for that
ps= merge_phyloseq(ps, repseqs)

# Metadata
#importing the metadata
metadata <- read.table("./metadata.tsv", header = TRUE, sep = "\t")

#converting non numeric and non-logical variables to factors
for(i in seq_len(ncol(metadata))) {
 if(!is.numeric(metadata[[i]]) && !is.logical(metadata[[i]]) && !is.integer(metadata[[i]])) {
    metadata[[i]] = as.factor(metadata[[i]])  } else {
     metadata[[i]]
 } 
}

#changing sample names
asvs <- otu_table(ps) %>% as.matrix

#Merging all artifacts
pst = phyloseq(otu_table(asvs, taxa_are_rows = TRUE), phy_tree(ps), sample_data(metadata), refseq(ps), tax_table(ps)) 
```

将 Illumina 输出导入 r

```r
ft = read_tsv('./MAC2023.github.io/data/illumina/Results_MAC23_MAC96_S96/OTU-tables/zOTU_table_GG.txt' ) %>% data.frame()
ft <- column_to_rownames(ft, "X.OTU.ID")

tx = ft %>% as.data.frame() %>% select(81)

ft[,dim(ft)[2]] <-NULL

tx = tidyr::separate(tx, col = taxonomy, into =c("Kingdom", "Phylum", 
            "Class", "Order", "Family", 
            "Genus", "Species"), 
            sep = ';') %>% as.matrix()

tx = apply(tx, 2, function(x) {gsub("^.__", "", x)})

colnames(ft) <- sapply(colnames(ft), function(x) {gsub(".*[S]([0-9]+)$", "BRK\\1", x)}) %>% as.vector()

mt = read.table('./data/metadata.txt', header = TRUE) 
mt
data.frame(barcode = colnames(ft))

mt[!rownames(mt) %in% colnames(ft),] %>% select(barcode) %>% pull()

tr$tip.label %>% length()
tr = ggtree::read.tree('./MAC2023.github.io/data/illumina/Results_MAC23_MAC96_S96/zOTU.tree') 

phyloseq(otu_table(ft, TRUE), phy_tree(tr), tax_table(tx)) %>% ggtree(layout = "circular", aes(color = Phylum))
```

将 ONT 输出导入 R

```r
#setwd('/home/hackion/Dropbox/Old-2gb/postdoc/mac2023')

#Reading feature table
ft = read_tsv('./data/count_matrix.tsv') %>% 
as.data.frame() %>% 
column_to_rownames('#OTU ID') %>% 
as.matrix()


# Reading taxonomy table
tx = read.table('./data/taxonomy.tsv', 
                header = F, 
                sep = '\t', 
                row.names = 1)
dim(tx)
```

```r
## [1] 5056    1
```
```r
# Parse the taxonomy
tx = tidyr::separate(data = tx, 
            col = V2, 
            into = c("Kingdom", "Phylum", 
            "Class", "Order", "Family", 
            "Genus", "Species"), 
            sep = ';') 

#removing taxon tags 

tx = apply(tx, 2, function(x) {gsub("^.__", "", x)})


#reading metadata
mt = read.table('./data/metadata.txt', header = TRUE) 


# Reading rep_seqs
seqs = Biostrings::readDNAStringSet('./data/rep_seqs.fasta')

# Reading the tree
tr = ggtree::read.tree('./data/tree.nwk')


# Merging data
pst = phyloseq(otu_table(ft, taxa_are_row=TRUE), phyloseq::tax_table(tx), phy_tree(tr), sample_data(mt), refseq(seqs))

pst
```

```r
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5056 taxa and 95 samples ]
## sample_data() Sample Data:       [ 95 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 5056 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 5056 tips and 4980 internal nodes ]
## refseq()      DNAStringSet:      [ 5056 reference sequences ]
```
```r
sum(rowSums(pst@otu_table)==0)
```
```r
## [1] 0
```

#### 3、reads的过滤和预处理

```r
#removing unassigned and NA phylum taxa
pst <- subset_taxa(pst, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "unassigned"))

#keeping only bacterial kingdom
pst <- subset_taxa(pst, Kingdom %in% "Bacteria")

pst
```
```r
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4858 taxa and 95 samples ]
## sample_data() Sample Data:       [ 95 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 4858 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 4858 tips and 4784 internal nodes ]
## refseq()      DNAStringSet:      [ 4858 reference sequences ]
```

基于（低）流行率的分类过滤：监督

```r
#monitoring the number of the samples in which the prevalence of a taxon is at least one
prevdf <- apply(otu_table(pst),ifelse(taxa_are_rows(pst), 1, 2), function(x){sum(x>0)})

weight_df <- data.frame(ASVprev = prevdf, 
                    TaxaAbund = taxa_sums(pst),
                    tax_table(pst))
head(weight_df)
```
```r
##          ASVprev TaxaAbund  Kingdom           Phylum          Class
## OTU_1133      14        22 Bacteria Actinobacteriota Actinobacteria
## OTU_3649      28        50 Bacteria Actinobacteriota Actinobacteria
## OTU_3778      15        30 Bacteria Actinobacteriota Actinobacteria
## OTU_614       70       934 Bacteria Actinobacteriota Actinobacteria
## OTU_1186      51       110 Bacteria       Firmicutes     Clostridia
## OTU_4237      22        35 Bacteria       Firmicutes     Clostridia
##                      Order             Family           Genus
## OTU_1133 Bifidobacteriales Bifidobacteriaceae Bifidobacterium
## OTU_3649 Bifidobacteriales Bifidobacteriaceae Bifidobacterium
## OTU_3778 Bifidobacteriales Bifidobacteriaceae Bifidobacterium
## OTU_614  Bifidobacteriales Bifidobacteriaceae Bifidobacterium
## OTU_1186    Lachnospirales    Lachnospiraceae         Blautia
## OTU_4237    Lachnospirales    Lachnospiraceae         Blautia
##                               Species
## OTU_1133 Bifidobacterium_adolescentis
## OTU_3649       Bifidobacterium_longum
## OTU_3778      Bifidobacterium_bifidum
## OTU_614  Bifidobacterium_adolescentis
## OTU_1186                         <NA>
## OTU_4237                         <NA>
```

```r
#Find out the phyla that are of mostly low-prevalence features by computing the total and average prev of features in Phylum
plyr::ddply(weight_df, "Phylum", function(x){cbind(means = round(mean(x$ASVprev), 2), sums = round(sum(x$ASVprev),2))}) %>% mutate(sanity = ifelse(means == sums, "TRUE", "FALSE"))
```

```r
##              Phylum means   sums sanity
## 1  Actinobacteriota 40.19  17645  FALSE
## 2      Bacteroidota 29.43  25282  FALSE
## 3     Cyanobacteria 20.25    486  FALSE
## 4  Desulfobacterota 22.89    435  FALSE
## 5        Firmicutes 35.80 116352  FALSE
## 6    Fusobacteriota 24.02    985  FALSE
## 7    Proteobacteria 26.01   4941  FALSE
## 8     Spirochaetota  9.80     49  FALSE
## 9 Verrucomicrobiota 35.52   1101  FALSE
```

根据所有样本中 n 个样本中出现的流行阈值过滤 ASV

```r
# This means that each ASV should have appeared at least in n samples to be kept.
asv.filter = function(asvtab, n.samples = 1 ){
  filter.threshold <- n.samples/ncol(asvtab) *100 # In how many samples out of total samples an ASV should have occured
  table_count <- apply(otu_table(asvtab), 2, function(x) ifelse(x>0, 1, 0)) %>% as.data.frame()
  suspected_ASV = table_count[which((rowSums(table_count)/ncol(table_count))*100 < filter.threshold),] %>% rownames()
  
  return(suspected_ASV)
}

print(glue("{ asv.filter(asvtab = otu_table(pst), n.samples = 2) } occured only in two samples"))

sus_ASV = asv.filter(asvtab = otu_table(pst), n.samples = 2)      
## Will you remove them?
#pst = subset_taxa(pst, !taxa_names(pst) %in% sus_ASV)


#or you could also do it by the phyloseq function

condition <- function(x) x>0
TaxaTokeep <- genefilter_sample(pst,condition,2)                      
#pst = subset_taxa(pst, taxa_names(pst) %in% TaxaTokeep)
```

根据丰度删除单调
```r
#A function to find singletones. You need to be careful about this step!
out.ASV = function(phyloseq, threshold =1, binwidth = 0.01) {
  
#Loading necessary pkgs      
  pacman::p_load(glue, tidyverse, reshape2, ggrepel, S4Vectors) # nolint
#This function requires phyloseq, tidyverse and glue packages to be loaded. 
    if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 100 ) {#making the relative abundance table
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
    } else if (sum(colSums(otu_table(phyloseq)))/ncol(otu_table(phyloseq)) == 1) {
                    rel_abund = as(t(otu_table(phyloseq)), "matrix")
                    } else {
                    rel_abund = as(t(apply(otu_table(phyloseq), 
                    ifelse(taxa_are_rows(phyloseq), 1,2), 
                    function(x) x/sum(x))), "matrix")  
                    } 
                      
                      
      names.single = apply(rel_abund, 1, function(x){ifelse(x == threshold, TRUE, 
                                                    ifelse(x == sum(x), TRUE, FALSE))}) %>% reshape2::melt() %>% 
                                                    filter(value == TRUE) %>% 
                                                    dplyr::select(2) %>% 
                                                    pull%>% as.vector()
                      
    if (length(names.single) == 0 ) {
          print(glue("WOW! {length(names.single)} singletones detected in this dataset"))
                 qplot.noSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
                  show.legend = F, main = "Frequency count of relative abundance, no singletones detected") +
         xlab ("Relative abundance in samples") + 
                   ylab("Frequency") + theme_bw()
      
    return(structure(list(qplot.noSing)))
                            
                        } else { 
                             
single.ASV = rel_abund[rownames(rel_abund) %in% names.single,]
single.ASV[single.ASV == 0] <- NA # A separate dataset for annotation of singletones on the barplot
                            
         qplot.withSing = qplot(rel_abund, geom = "histogram", binwidth = binwidth, 
        main = "Frequency count of relative abundance with singletones") +
        geom_bar(aes(single.ASV), fill = "red",  color = NA, width = binwidth)+
                       xlab ("Relative abundance in samples") + ylab("Frequency") + 
                       geom_label_repel(aes(x = 1, y =length(rel_abund)/5), 
                       label.padding =  unit(0.55, "lines"), 
                       label = glue("{length(names.single)}\n Singletones"), color = "black") + 
           theme_bw()
                            
                       qplot.rmSing = qplot(rel_abund[!rownames(rel_abund) %in% names.single, ], geom = "histogram",
                       binwidth = binwidth, main = "Frequency count of relative abundance without singletones") +
                       xlab ("Relative abundance in samples") + ylab("Frequency")+ theme_bw()
                            
                       print(glue('Oh no..! {length(names.single)} singletones detected in the dataset'))
                       return(structure(list(qplot.withSing, qplot.rmSing, unlist(names.single))) )
                    
                        }                        
    
                             
        }
                        
single.test = out.ASV(phyloseq = pst, threshold = 2, binwidth = 0.01)
#singletones = single.test[[3]] #here you can extract the names of the singletones

single.test[[1]]#to show the plot with singletones
single.test[[2]]#to show the plot without singletones

#Now you can remove the singletones from your pst file as follows:
#pst = subset_taxa(pst, !taxa_names(ps)%in% singletones)


rm(single.test)
```

#### 4.阿尔法多样性

##### 4.1. 稀疏

```r
#library(MicrobiotaProcess)


#This takes a bit of time
ps_rar_curve <- MicrobiotaProcess::ggrarecurve(obj = pst, 
                    indexNames = c("Observe", "Shannon"),
                    chunks=400, 
                    theme(legend.spacing.y = unit(0.02, "cm"),
                          legend.text = element_text(size = 6)), 
                          show.legend=F) 



ps_rar_curve + theme_bw() + geom_vline(xintercept = 30000, lty = 2, color = alpha("red", 0.5)) + 
  ggtitle("Rarefaction curves based on Richnessh and evenness", "Each line is a sample")

```

```r
#ggsave("./Alpha/rarefaction.curve.jpeg", device = "png", dpi = 300, height = 6, width = 9)

#rarefying the table with minimum depth of 10000 reads per sample
ps_rar = rarefy_even_depth(pst, sample.size = 30000, replace = FALSE)#in this depth we have lost one sapmle (F.29) and no ASVs.

# Taxonomic filtering based on abundance for rarefied data: supervsed
## Abumdance: ASV > 0.0001% overall abundance across all samples
total.depth <- sum(otu_table(ps_rar))
totAbuThreshold <- 1e-4 * total.depth
ps_rar <- prune_taxa(taxa_sums(ps_rar)>totAbuThreshold, ps_rar)

ps_rar
```

```r
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1221 taxa and 87 samples ]
## sample_data() Sample Data:       [ 87 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 1221 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1221 tips and 1214 internal nodes ]
## refseq()      DNAStringSet:      [ 1221 reference sequences ]
```

```r
sum(rowSums(ps_rar@otu_table)==0)
```

```r
## [1] 0
```

##### 4.2. Alpha 多样性指标

```r
#Calculating the alpha diversity indexes

## Richness
Chao1 = estimate_richness(pst, split = TRUE, measures = "Chao1")#for richness, we don't use rarefied table



## Evenness: shannon => H' = -Σ (Pi * ln(Pi))
Shannon = estimate_richness(ps_rar, split = TRUE, measures = "Shannon")

#Faith Phylogenetic Diversity: PD = Σ(branch lengths leading to selected species)
library(picante, verbose = FALSE)
FaithPD = picante::pd(samp = t(otu_table(ps_rar)), tree = phy_tree(ps_rar))$PD

#adding the indexes to the metadatas              
sample_data(ps_rar) <- data.frame(sample_data(ps_rar), Chao1=Chao1[rownames(Chao1)%in% rownames(sample_data(ps_rar)),][[1]], Shannon = Shannon$Shannon, FaithPD = FaithPD) #note that we have removed that sample which has been removed by rarefaction.  
```

##### 4.3.可视化α多样性指数

```r
library(ggpubr, verbose = FALSE)
library(reshape2, verbose = FALSE)

alpha_ccd = sample_data(ps_rar)   %>% data.frame


long_mtdat <- melt(alpha_ccd, )
long_mtdat<- long_mtdat[long_mtdat$variable %in% c("Chao1", "Shannon", "FaithPD"),]



long_mtdat$variable <- factor(long_mtdat$variable , levels =  c("Chao1", "Shannon", "FaithPD"))

pie(rep(10, 9), col = colors()[c(190, 20, 30, 15, 8, 90, 35, 120, 42)], clockwise = TRUE, labels = as.character(colors()[c(1, 20, 30, 15, 8, 90, 35, 120, 42)]))
```
```r
cols = colors()[c(190, 20, 30, 15, 8, 90, 35, 120, 42)]

my_comp = list(c("A", "B"), c("A", "C"), c("A", "E"), c("A", "Mock"), c("A", "B"))
#comparison between Samples by SampleID and digesta

alpha_p = ggplot(long_mtdat, aes(x = group, y = value)) +
  geom_violin(aes(fill = group), trim = F) +
  stat_compare_means(paired = FALSE, comparison = my_comp, method = "t.test", label = "p.signif") +
  geom_boxplot(width = 0.10) + 
  geom_jitter(color = "black", alpha = 0.5)+ 
  facet_wrap(~long_mtdat$variable, scales = "free_y") + 
  theme_bw() + 
  scale_fill_manual(values = cols) + 
  theme(legend.title = element_text( size = 15, face = "bold"), 
  axis.title.x = element_text( face = "bold", size = 15), axis.text.x = element_text( size = 15), 
  axis.title.y = element_text( face = "bold", size = 15),
  axis.text.y = element_text( size = 15), strip.text.x = element_text( size = 15, face = "bold")) + 
  labs(fill = "Sample type",  y = "Alpha diversity",
  title = "Alpha diversity metrics in different group samples.\n Means are compared by unpaired t.test") + 
  xlab("Sample type")


#ggsave(plot = alpha_p, "./Alpha/alpha.jpeg", device = "jpeg", width = 15, dpi =300)

alpha_p
```
```r
lm1 = lm(Chao1 ~ group, data = alpha_ccd)

summary(lm1)
```
```r
## 
## Call:
## lm(formula = Chao1 ~ group, data = alpha_ccd)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -679.86 -196.49   40.58  234.99  533.47 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  2588.961     92.412  28.015  < 2e-16 ***
## groupB        335.459    130.690   2.567 0.012177 *  
## groupC        -31.262    133.628  -0.234 0.815637    
## groupD          5.494    130.690   0.042 0.966577    
## groupE        535.158    130.690   4.095 0.000102 ***
## groupF        169.873    130.690   1.300 0.197493    
## groupG       -802.872    170.400  -4.712 1.05e-05 ***
## groupH      -1066.241    141.162  -7.553 6.86e-11 ***
## groupMock   -1103.125    244.499  -4.512 2.24e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 320.1 on 78 degrees of freedom
## Multiple R-squared:  0.7232, Adjusted R-squared:  0.6948 
## F-statistic: 25.47 on 8 and 78 DF,  p-value: < 2.2e-16
```

#### 5. 贝塔多样性

计数数据的转换又称为特征缩放

```r
library("patchwork", verbose = FALSE)

#Untransformed
p1 = ggplot() +  
  geom_histogram(aes(x = rowSums(otu_table(pst))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 2000) +
  theme_bw() + 
    ggtitle("Histogram for raw") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p2 = ggplot() +  
  geom_histogram(aes(x = log10(rowSums(otu_table(pst)))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 0.05) +
  theme_bw() + 
    ggtitle("Histogram for raw, Log10T") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p1 + p2 + plot_layout(ncol = 2 , heights = 8)
```
```r
#Rarefied
p3 = ggplot() +  
  geom_histogram(aes(x = rowSums(otu_table(ps_rar))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 2000) +
  theme_bw() + 
    ggtitle("Histogram for rarefied data") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p4 = ggplot() +  
  geom_histogram(aes(x = log10(rowSums(otu_table(ps_rar)))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 0.07) +
  theme_bw() + 
    ggtitle("Histogram for raw, rarefied Log10T") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p3 + p4 +plot_layout(ncol = 2 , heights = 8)
```
```r
#Relabund
ps_rel = transform_sample_counts(ps_rar, function(x) {x / sum(x) * 100})

p5 = ggplot() +  
  geom_histogram(aes(x = rowSums(otu_table(ps_rel))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 5) +
  theme_bw() + 
    ggtitle("Histogram for relabund data") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p6 = ggplot() +  
  geom_histogram(aes(x = log10(rowSums(otu_table(ps_rel)))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 0.07) +
  theme_bw() + 
    ggtitle("Histogram for relabund, Log10T") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

 p5 + p6+ plot_layout(ncol = 2 , heights = 8)
```
```r
#Natural log
ps_log = transform_sample_counts(ps_rar, function(x) (log(1+x)))

p7 = ggplot() +  
  geom_histogram(aes(x = rowSums(otu_table(ps_log))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 5) +
  theme_bw() + 
    ggtitle("Histogram for log data") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p8 = ggplot() +  
  geom_histogram(aes(x = log10(rowSums(otu_table(pst)))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 0.07) +
  theme_bw() + 
    ggtitle("Histogram for log, Log10T") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p7 + p8 + plot_layout(ncol = 2 , heights = 8)
```
```r
#VST
vst_count =  varianceStabilizingTransformation(object = as(otu_table(ps_rar),  "matrix"), fitType = "mean", blind = TRUE)

ps_vst = ps_rar
otu_table(ps_vst) <- otu_table(vst_count, taxa_are_rows = TRUE)

p9 = ggplot() +  
  geom_histogram(aes(x = rowSums(otu_table(ps_vst))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 10) +
  theme_bw() + 
    ggtitle("Histogram for VST data") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p10 = ggplot() +  
  geom_histogram(aes(x = log10(rowSums(otu_table(ps_vst)))), 
  fill = "#4a8d78", color = "#ffffff", binwidth = 0.07) +
  theme_bw() + 
    ggtitle("Histogram for VST, Log10T") + 
    xlab("Taxa sum") + 
    ylab("Count, all samples")

p9 + p10 + plot_layout(ncol = 2 , heights = 8)
```
```r
msd = meanSdPlot(otu_table(ps_rar), plot = FALSE, ranks = TRUE)

p11 = msd$gg + ggtitle("Mean vs. STD for raw counts of taxa") +
theme_bw() + 
xlab("Mean of normalized raw counts") + 
ylab("Standard deviation of counts") +
guides(fill = "none")


msd2 = meanSdPlot(otu_table(ps_vst), plot = FALSE, ranks = TRUE)

p12 = msd2$gg + 
ggtitle("Mean vs. STD for VST counts of taxa") +
theme_bw() + 
xlab("Mean of normalized VST counts") + 
ylab("Standard deviation of VST counts")

p11 + p12 + plot_layout(ncol = 2 , heights = 8)
```
```r
# Standardization  (x-mean(x)/sd(x))
tb = ps_rar@otu_table %>% data.frame()
std = apply(tb, 2, function(x){scale(x)})

rownames(std) <- rownames(tb)

qplot(rowSums(std), bins = 100) 
```
```r
pca_result = prcomp(std, scale = FALSE)
names(pca_result)
```
```r
## [1] "sdev"     "rotation" "center"   "scale"    "x"

pca_result$center

##         BRK01         BRK02         BRK03         BRK04         BRK05 
## -4.841882e-18 -4.516532e-18 -1.260552e-17  3.100055e-18 -3.685399e-18 
##         BRK06         BRK07         BRK08         BRK09         BRK10 
##  7.618008e-18  2.040467e-17  3.316007e-18 -1.672069e-17 -6.933211e-19 
##         BRK11         BRK12         BRK13         BRK14         BRK15 
## -4.578867e-18  1.868983e-17 -9.308688e-18 -9.223444e-18  6.710865e-18 
##         BRK16         BRK17         BRK18         BRK19         BRK20 
##  1.265311e-17  3.233746e-17  2.415258e-19  8.945689e-18  1.014977e-17 
##         BRK21         BRK22         BRK23         BRK24         BRK25 
## -3.020138e-18 -2.631637e-17 -1.553153e-17 -1.329742e-17 -1.175662e-18 
##         BRK26         BRK27         BRK28         BRK29         BRK30 
##  7.176868e-18 -4.081786e-18 -7.600248e-18 -1.728401e-17  1.236790e-17 
##         BRK31         BRK32         BRK33         BRK34         BRK36 
##  1.065555e-17 -5.702317e-18  1.599895e-17 -9.921737e-18 -1.142879e-17 
##         BRK37         BRK38         BRK39         BRK40         BRK41 
## -1.132330e-18  4.797839e-18 -1.785160e-17 -9.740593e-18 -1.431395e-18 
##         BRK42         BRK43         BRK44         BRK45         BRK46 
##  7.954723e-18  2.809868e-18  9.705075e-18  1.118122e-18 -1.449013e-17 
##         BRK47         BRK48         BRK49         BRK50         BRK51 
## -1.050495e-17  3.361471e-18  8.628154e-18 -6.707313e-18  6.816000e-18 
##         BRK52         BRK53         BRK54         BRK55         BRK56 
##  1.298627e-17 -3.862992e-18 -2.711482e-18  4.148383e-18  1.683044e-18 
##         BRK57         BRK58         BRK59         BRK60         BRK61 
##  1.528432e-17 -1.047796e-17 -2.039188e-17 -1.482258e-17 -1.254087e-17 
##         BRK62         BRK63         BRK64         BRK65         BRK66 
##  3.474420e-18  3.078743e-17 -1.463362e-17 -1.260125e-17 -1.770952e-18 
##         BRK67         BRK68         BRK69         BRK70         BRK71 
## -2.333565e-18  6.753132e-18  1.052697e-17  8.008267e-18  4.469648e-18 
##         BRK72         BRK81         BRK82         BRK83         BRK84 
##  1.414045e-17  1.019807e-17 -3.789824e-18 -3.597314e-18 -9.277432e-18 
##         BRK85         BRK86         BRK87         BRK88         BRK89 
##  9.331420e-18  4.660027e-19  1.505984e-17 -7.615166e-18 -1.365615e-17 
##         BRK90         BRK91         BRK92         BRK93         BRK94 
##  1.325266e-17  9.757642e-18 -3.472288e-18 -8.391956e-18 -1.005031e-17 
##         BRK95         BRK96 
##  6.060166e-18  1.130554e-17
```
```r
pca_df = data.frame(PC1 = pca_result$rotation[,1], PC2 = pca_result$rotation[,2], group = sample_data(ps_rar)$group)


explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc_labels <- paste0("PC", 1:length(explained_variance))
var_df <- data.frame(PC = pc_labels, ExplainedVariance = explained_variance)


ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
geom_point(size = 8, alpha = 0.6) +
scale_color_manual(values = colors()[c(190, 40, 30, 15, 8, 90, 35, 80, 42)]) + theme_bw() +
xlab(glue("PC1, {round(var_df$ExplainedVariance[1] * 100,1)}%")) +
ylab(glue("PC1, {round(var_df$ExplainedVariance[2] * 100,1)}%")) +
guides(color = guide_legend(title = "Group")) +
geom_vline(xintercept = 0, lty = 2, color = alpha("orange", alpha = 0.7))+
geom_hline(yintercept = 0, lty = 2, color = alpha("orange", alpha = 0.7))
```

##### 5.1. 排序图

```r
cols = colors()[c(190, 20, 30, 15, 8, 90, 35, 120, 42)]
#bray 

#bray PCoA: BC = (Σ|X_i - Y_i|) / (Σ(X_i + Y_i))
bray_pcoa=ordinate(ps_log, method="PCoA", distance = "bray")
evals<-bray_pcoa$values$Eigenvalues

bray_p = plot_ordination(ps_log, bray_pcoa,  title = "Bray-Curtis PCoA plot, Log", color = "white") +
  labs(col="Group")+ geom_point( aes(fill = group), alpha = 0.75, pch = 21, size = 8,  color = "black", show.legend = FALSE) +
   #coord_fixed() + 
  #stat_ellipse(aes(group = group, fill = group), 
  #  show.legend = F, type = "t", level = 0.05, lty = 2, geom = "polygon", alpha = 0.1) + 
    labs(x = sprintf("PCo1 [%s%%]", round(evals/sum(evals)*100,1)[1]),
      y = sprintf("PCo2 [%s%%]", round(evals/sum(evals)*100, 2)[2])) +
         scale_fill_manual(values = cols) +
            geom_vline(xintercept = 0, lty = 2, alpha = 0.5, color = "blue") +
             geom_hline(yintercept = 0, lty = 2, alpha = 0.5, color = "blue") + theme_bw()+ 
              theme(axis.title = element_text(face = "bold"), 
                legend.title = element_text(size = 10, face = "bold"),
                 legend.text = element_text(face = "bold"), 
                   axis.text = element_text(size = 15))


#ggsave("./Beta/bray.pcoa.dig.vs.muc.jpeg", dpi = 300)

wunifrac = ordinate(ps_log, method="NMDS", distance='wunifrac')
```
```r
## Run 0 stress 0.04758852 
## Run 1 stress 0.05257593 
## Run 2 stress 0.04758856 
## ... Procrustes: rmse 3.155827e-05  max resid 8.339069e-05 
## ... Similar to previous best
## Run 3 stress 0.04758858 
## ... Procrustes: rmse 3.471801e-05  max resid 8.939346e-05 
## ... Similar to previous best
## Run 4 stress 0.0475886 
## ... Procrustes: rmse 1.71134e-05  max resid 5.518678e-05 
## ... Similar to previous best
## Run 5 stress 0.04758853 
## ... Procrustes: rmse 1.978362e-05  max resid 5.187433e-05 
## ... Similar to previous best
## Run 6 stress 0.07429898 
## Run 7 stress 0.05257595 
## Run 8 stress 0.07429908 
## Run 9 stress 0.04758852 
## ... Procrustes: rmse 1.726039e-05  max resid 4.740674e-05 
## ... Similar to previous best
## Run 10 stress 0.04758851 
## ... New best solution
## ... Procrustes: rmse 1.241057e-05  max resid 2.7403e-05 
## ... Similar to previous best
## Run 11 stress 0.04758852 
## ... Procrustes: rmse 4.215754e-06  max resid 1.148943e-05 
## ... Similar to previous best
## Run 12 stress 0.07429939 
## Run 13 stress 0.04758853 
## ... Procrustes: rmse 8.518784e-06  max resid 2.341394e-05 
## ... Similar to previous best
## Run 14 stress 0.04758852 
## ... Procrustes: rmse 7.077508e-06  max resid 1.774362e-05 
## ... Similar to previous best
## Run 15 stress 0.04758853 
## ... Procrustes: rmse 1.264987e-05  max resid 3.6171e-05 
## ... Similar to previous best
## Run 16 stress 0.04758853 
## ... Procrustes: rmse 3.895438e-05  max resid 0.0001050886 
## ... Similar to previous best
## Run 17 stress 0.04758853 
## ... Procrustes: rmse 4.289572e-05  max resid 0.0001150774 
## ... Similar to previous best
## Run 18 stress 0.05257595 
## Run 19 stress 0.04758852 
## ... Procrustes: rmse 1.454469e-06  max resid 7.042092e-06 
## ... Similar to previous best
## Run 20 stress 0.04758853 
## ... Procrustes: rmse 6.52229e-06  max resid 2.010471e-05 
## ... Similar to previous best
## *** Best solution repeated 9 times
```
```r
wunifrac_p = plot_ordination(ps_log, wunifrac, title = "WUNIFRAC NMDS, LogT") +
   geom_point( aes(fill = group), pch = 21, 
   alpha = 0.75, color = "black", size =8, show.legend = TRUE)  + 
   #coord_fixed()+
    labs(fill="Group")+
         scale_fill_manual(values = cols) +
            geom_vline(xintercept = 0, lty = 2, alpha = 0.5, color = "blue") +
             geom_hline(yintercept = 0, lty = 2, alpha = 0.5, color = "blue") + theme_bw()+ 
              theme(axis.title = element_text(face = "bold"), 
                legend.title = element_text(size = 10, face = "bold"),
                 legend.text = element_text(face = "bold"), 
                   axis.text = element_text(size = 15))

bray_p + wunifrac_p + plot_layout( heights = 8, widths = 15, ncol = 2) + plot_annotation("Ordination plots of compostional analysis", "Bray-Curtis dissimilarity and\n Weighted Unifrac phylogenetic distance", "You might know this as Beta diversity")
```

#### Gloomer：一个包装“tax_glom()”函数并向分类单元添加简洁唯一名称的函数

```r
# Gloomer

#A function to create unique names for each ASV. If species is set as the taxa level, it removes any NA in Order level then attempts to use the name of one level higher taxa for those who have similar names, e.g. uncultured_bacterium

gloomer = function(ps = data, taxa_level = taxa_level, NArm = "TRUE"){
    rank.names = c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    

#====================Sometimes in genus level, we might have multiple uncultured organisms, which if we want to make unique out of them for the species level it won't work====
    #since adding uncultured to uncultered is sill duplication. therefore if the taxa_level is set to species we first make a unique genus and then we go further to the speices===#

#Removing unculured Family
ps = subset_taxa(ps, !Family %in% c("uncultured", "NA", "uncategorized", "unassigend", "", " "))
    
if(taxa_level == "Species") {

    ps = subset_taxa(ps, !Genus %in% NA)#we remove genus tagged NA
tax_table(ps)[, taxa_level] <- ifelse(is.na(tax_table(ps)[, taxa_level]), paste0("unknown"), paste(tax_table(ps)[, taxa_level]))#convert NA in species into unknown
    
  physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
  taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]

   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
   
#first take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] %in% c("uncategorized", NA, "uncultured", "unassigend", "", " "),
       paste0("[", taxdat[,length(rank.names[1:which(rank.names=="Genus")])-1], "]", "_", taxdat[,6]), taxdat[,6])
    
spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") 

if(dim(duplis)[[1]] > 0) {
duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1) %>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis,
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level]))

#check if all the names are unique at species level, otherwise we will bring family instead of genus
   dupies <-  taxdat[duplicated(taxdat[,"uni"]), "uni"] 
    if(length(dupies)>0) {
        taxdat = taxdat %>% data.frame %>% mutate( uni2= ifelse(taxdat[, "uni"] %in% dupies,
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-2], "]", "_", taxdat[,"uni"]), taxdat[,"uni"]))
        
        taxdat[, taxa_level] = taxdat[, "uni2"]
        taxdat[, "uni"] <- NULL
        taxdat[, "uni2"] <- NULL
        taxdat <- as(taxdat, "matrix")   
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[, taxa_level]
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
        
    }
    else 
    {
        
taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as(taxdat, "matrix")   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
           }
    
} else {
    
taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
    
}
       
    
#==========================================# 
} else if (taxa_level == "Genus") {
    
    physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = NArm)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
   taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
    otudat = otu_table(physeq)
    
# take care of the uncultured genus
taxdat[,6] = ifelse(taxdat[,6] %in% c("uncategorized", NA, "uncultured", "unassigend", "", " "),
       paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level])
    
gen1 = taxdat[, taxa_level] %>% as.vector
gen2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(gen2), nrow = length(gen1))
    for(i in seq_along(gen1)){
        for(j in seq_along(gen2)){
    uni[i, j] = ifelse(gen1[i] == gen2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-gen1
colnames(uni) <- gen2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

        if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name
    
        duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
        taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste0("[", taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "]", "_", taxdat[,taxa_level]), taxdat[,taxa_level]))
    
        taxdat[, taxa_level] = taxdat[, "uni"]
        taxdat[, "uni"] <- NULL

        taxdat <- as(taxdat, "matrix")
 
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[taxdat[,taxa_level] %in% rownames(otudat), taxa_level]
        taxdat <- as.matrix(taxdat) 
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
 
        } else {

        taxdat <- as.matrix(taxdat) 
        taxdat <- tax_table(taxdat)
        rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
        rownames(taxdat) <- taxdat[, taxa_level]
        taxdat <- tax_table(taxdat)
        taxa_names(physeq) <- taxa_names(taxdat)
        tax_table(physeq) <- taxdat
        otu_table(physeq) <- otudat
       }   
    
} else {
    
    
physeq = tax_glom(physeq = ps, taxrank = taxa_level, NArm = TRUE)
    taxdat = tax_table(physeq)[, seq_along(rank.names[1:which(rank.names == taxa_level)])]
    
taxdat = taxdat[complete.cases(taxdat),] %>% as.data.frame
otudat = otu_table(physeq)
    
spec1 = taxdat[, taxa_level] %>% as.vector
spec2  = taxdat[, taxa_level] %>% as.vector

    uni  = matrix(NA, ncol = length(spec2), nrow = length(spec1))
    for(i in seq_along(spec1)){
        for(j in seq_along(spec2)){
    uni[i, j] = ifelse(spec1[i] == spec2[j] , "TRUE", "FALSE")
    }
        }

rownames(uni) <-spec1
colnames(uni) <- spec2   
uni[upper.tri(uni, diag = TRUE)] = 0 #get rid of diagonals and upper triangle

duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE")

if(dim(duplis)[[1]] > 0){#if there is not duplications, we can simply use the taxa names as the row name
    
    duplis = uni %>% reshape2::melt() %>% filter(value == "TRUE") %>% dplyr::select(1)%>% unique() %>% unlist %>% as.vector
taxdat = taxdat %>% mutate( uni= ifelse(taxdat[, taxa_level] %in% duplis, 
                    paste(taxdat[,length(rank.names[1:which(rank.names==taxa_level)])-1], "_", taxdat[,taxa_level]), taxdat[,taxa_level]))
    
taxdat[, taxa_level] = taxdat[, "uni"]
taxdat[, "uni"] <- NULL
taxdat <- as.matrix(taxdat)   
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
} else {

taxdat <- as.matrix(taxdat) 
taxdat <- tax_table(taxdat)
rownames(otudat) <- taxdat[rownames(taxdat) %in% rownames(otudat), taxa_level]
rownames(taxdat) <- taxdat[, taxa_level]
taxdat <- tax_table(taxdat)
taxa_names(physeq) <- taxa_names(taxdat)
tax_table(physeq) <- taxdat
otu_table(physeq) <- otudat
}
#ps = phyloseq(otu_table(otudat, taxa_are_rows = T), tax_table(as.matrix(taxdat)), sample_data(physeq))
 

}
return(physeq) 
    }
```

##### 5.2. 基于图表的 beta 多样性分析

```r
library(igraph, verbose = FALSE)
library(ggnetwork, verbose = FALSE)
library(phyloseqGraphTest, verbose = FALSE)
#For total none-rarefied dataset
pst.spec <- gloomer(ps = pst, taxa_level = "Species", NArm = TRUE)
ps_total = prune_taxa(taxa_sums(pst.spec)>1000, pst.spec)#filtering the taxa based on total sum
sample_data(ps_total)$sampleID <- rownames(sample_data(ps_total))

net <- make_network(ps_total, directed = FALSE, max.dist = 0.35, distance = "bray", type = "samples")

sampledata <- sample_data(ps_total) %>% data.frame

sampledata$sampleID <- rownames(sampledata)

V(net)$id <- sampledata[names(V(net)), "group"] %>% as.vector
V(net)$sample <-  rownames(sampledata)[rownames(sampledata) %in% names(V(net))] %>% as.vector
V(net)$sample_type <- sampledata[names(V(net)), "group"] %>% as.vector

#graph permutational test
set.seed(2023)
graph.test <- graph_perm_test(ps_total, sampletype = "group", max.dist = 0.35,
              grouping = "sampleID", distance = "bray", type = "mst", 
                             nperm = 1000)
                      
print(glue("Test statistic for the graph-based analysis indicates that\n with p-value = {graph.test$pval}, \n we REJECT the Null hypothesis that the distribution of taxa in samples from different groiups are similar!"))
```
```r
## Test statistic for the graph-based analysis indicates that
## with p-value = 0.000999000999000999, 
## we REJECT the Null hypothesis that the distribution of taxa in samples from different groiups are similar!
```
```r
V(graph.test$net)$group <- sampledata[names(V(graph.test$net)), "group"] %>% as.vector

plotNet1 = phyloseqGraphTest::plot_test_network(graph.test) + 
theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9)) + 
geom_nodes(size = 4, aes(color = sampletype))+
scale_color_manual(values = cols)


plotPerm1=plot_permutations(graph.test, bins = 40) + 
                  geom_text(aes(label = glue("{ifelse(round(graph.test$pval, 2) == 0,
                  'P < 0.01', round(graph.test$pval, 2))}"), 
                  x = 65, y =40), color = "red") + 
                  theme_bw()+ 
                  ggtitle("Permutation test for pure edges", "Bray-Curtis")

grid.arrange(ncol = 2,  plotNet1, plotPerm1) + geom_col(inherit.aes = F, color = "red") 
```
```r
## NULL

#ggsave("./graph_perm_total.jpeg", plot = net.grid, device = "jpeg", width = 15, beight = 10, dpi = 300)
```

##### 5.3. Beta多样性指数的统计分析：基于距离的冗余分析（dbRDA）

```r
#Calculating Bray-Curtis dissimilarity coefficeints
bray_log = phyloseq::distance(ps_log, method = "bray")
```

##### 5.3.1. 检查质心周围方差的离散程度（方差同源检验）

```r
#test for the disperssion of the variance around the centroids
library(vegan)
library(permute, verbose = FALSE)
#set the age of the animal as the random variable
set.seed(1990)
h <- with(data = data.frame(sample_data(ps_log)), how(blocks = group, nperm = 999))

##total data

#Now we do a Homogeneity of dispersion test
set.seed(10)
bray.disp <- vegan::betadisper(bray_log, group = sample_data(ps_log)$group, 
                                 type = "centroid")#if the p-value is significant, it means that there is a significant difference in variance for any of the tested levels. 
perm.test = permutest(bray.disp, permutation =h, pairwise = T)

p.val.perm = perm.test$tab$`Pr(>F)`[[1]]

disp.centroid = bray.disp$centroids %>% as.data.frame
disp.vectors = bray.disp$vectors %>% as.data.frame
eig.vals = bray.disp$eig



#jpeg( "./Beta/dispersion of variance_bray_total_sample.type.jpeg",  quality = 100)

plot(bray.disp, col = cols, 
    bty = "n",  las = 1, 
     main = "Dispersion of variance around the centroids, \n bray | LogT dataset", 
     sub=NULL,  xlab = sprintf("PCo1 [%s%%]", round(eig.vals/sum(eig.vals)*100,1)[1]),
      ylab = sprintf("PCo2 [%s%%]", 
    round(eig.vals/sum(eig.vals)*100,1)[2])); text(glue("P = {p.val.perm}"), 
        x = -0.1, y = -0.29, cex = 1.5, col = "red")
```

##### 5.3.2. 数据库RDA模型

```r
#Whole dataset
set.seed(2023)
h <- with(data = data.frame(sample_data(ps_log)), how(blocks = group, nperm = 999))

bray.dbrda = dbrda(t(otu_table(ps_log)) ~ group,
                        dist = "bray", permutations=h, data = sample_data(ps_log)%>%data.frame)
#Sex not sig, so reduced the model

bray.dbrda
```
```r
## Call: dbrda(formula = t(otu_table(ps_log)) ~ group, data =
## sample_data(ps_log) %>% data.frame, distance = "bray", permutations =
## h)
## 
##               Inertia Proportion Rank RealDims
## Total         13.5609     1.0000              
## Constrained    8.3821     0.6181    8        8
## Unconstrained  5.1788     0.3819   78       62
## Inertia is squared Bray distance 
## 
## Eigenvalues for constrained axes:
## dbRDA1 dbRDA2 dbRDA3 dbRDA4 dbRDA5 dbRDA6 dbRDA7 dbRDA8 
##  4.931  1.071  0.857  0.555  0.470  0.343  0.138  0.017 
## 
## Eigenvalues for unconstrained axes:
##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
## 1.9737 1.1133 0.6623 0.4507 0.1619 0.1129 0.0996 0.0569 
## (Showing 8 of 78 unconstrained eigenvalues)
```
```r
permutest(x = bray.dbrda, by = "terms", permutations = h)
```
```r
## 
## Permutation test for dbrda under reduced model 
## 
## Blocks:  group 
## Permutation: free
## Number of permutations: 999
##  
## Model: dbrda(formula = t(otu_table(ps_log)) ~ group, data =
## sample_data(ps_log) %>% data.frame, distance = "bray", permutations =
## h)
## Permutation test for all constrained eigenvalues
##          Df Inertia      F Pr(>F)
## group     8  8.3821 15.781      1
## Residual 78  5.1788
```

##### 5.3.3. 绘制模型提取图

```r
# Make the plot out of the model, which is the variation explained only by the terms 

#digesta
score.site = vegan::scores(bray.dbrda, display = "sites") %>% as.data.frame
score.centroid = vegan::scores(bray.dbrda, display = "cn")%>% as.data.frame
rownames(score.centroid)<- levels(sample_data(ps_log)$group)


eig.vals = bray.dbrda$CCA$eig
inertia.total = bray.dbrda$tot.chi #total variation (inertia) explained. 
#this number should be used as the denominator for measuring the amount of variance out of totoal variance wxplained by each dbrda

#Digesta
score.site %>% ggplot(aes(dbRDA1, dbRDA2)) +
geom_point(aes(fill = sample_data(ps_log)$group), color = "black", pch = 21, alpha = 0.5, size =6  ) + 
geom_hline(yintercept = 0, lty = 2, alpha =0.5) + 
geom_vline(xintercept = 0, lty = 2, alpha = 0.5)+
scale_fill_manual(values = cols) + 
theme_bw() + 
scale_y_continuous(na.value = c(-2, 3), n.breaks = 10) +
scale_x_continuous(na.value = c(-1, 1), n.breaks = 10) + 
labs(fill ="Groups") + 
xlab(label = paste("dbRDA1 [", round(eig.vals[[1]]/sum(eig.vals)*100, 1), 
                    "% of fitted and", 
                   round(eig.vals[[1]]/inertia.total*100, 1), 
                   "% of total variation]")) + 
ylab(label = paste("dbRDA2 [", round(eig.vals[[2]]/sum(eig.vals)*100, 1), 
                   "% of fitted and", round(eig.vals[[2]]/inertia.total*100, 1), 
                   "of total variation]")) + 
theme(axis.title = element_text(size = 10), 
      text = element_text(size = 13, face = "bold"),
      axis.text.x =element_text(size =10, face = "bold"),
axis.text.y =element_text(size =10, face = "bold")) + 
ggtitle(label = "dbRDA plot of Bray", "LogT data") 
```
```r
#ggsave("./Beta/bray.dbRDA.dig.jpeg", height = 8, width = 9, dpi =300)
```

#### 6. 差异丰度分析：DESeq2

##### 6.1. 堆积条形图

```r
#Aglomerating the taxa

phyl =  gloomer(ps_rar, taxa_level =  "Phylum", NArm = TRUE)
spec =  gloomer(ps_rar, taxa_level =  "Species", NArm = TRUE)

#Barplot of relative abundance of phylum in digesta
trans_phyl <- merge_samples(phyl, "group")
relabund_phyl <- transform_sample_counts(trans_phyl, function(x) x / sum(x)*100)
relabund_phyl = prune_taxa( taxa_sums(relabund_phyl)>0.0001, relabund_phyl)

#choosing colors
phylcol=c( "deepskyblue",'springgreen3','snow3','burlywood4', 'cadetblue', 'darkblue',
          'cornflowerblue','deeppink2','orangered',  'dimgrey', 'red','limegreen',
        'cyan1','darkmagenta', 'purple', 'cyan4', 'gold') 
phy_col = phylcol[1:7]

    
 phyl_p <- plot_bar(relabund_phyl, fill="Phylum") + 
    scale_fill_manual(values = phylcol) + 
    xlab("Groups") + 
    ylab("Relative abundance, %")+
    ggtitle(label = "Stacked barplot", "rarefied, relabund data > 0.01%")  +
            theme_bw() + 
            theme( legend.position = "right", 
            text = element_text(size =15, face = "bold")) #+ 
            #coord_flip()
                                             
 phyl_p                                            
```

##### 6.2. DESeq2

请记住使用原始数据，因为 DESEq2 内部会计算大小因子和色散的差异。

```r
phyl =  gloomer(pst, taxa_level =  "Phylum", NArm = TRUE)
spec = gloomer(pst, taxa_level = "Species", NArm = TRUE)

#converting phylosq to deseq
phyl_dds <- phyloseq_to_deseq2(phyl, design = ~  group) 
spec_dds <- phyloseq_to_deseq2(spec, design = ~ group)

#calculate geometric means prior to estimate size factors
gm.mean = function(x, na.rm= TRUE) {
    exp(sum(log(x[x>0]), na.rm=na.rm)/length(x))
}

##Phylum level
geo.mean = apply(counts(phyl_dds), 1, gm.mean)
phyl_dds = estimateSizeFactors(phyl_dds, geoMeans = geo.mean)
phyl_dds <-DESeq(phyl_dds, test = "Wald", fitType = "parametric")

#Species level
geo.mean = apply(counts(spec_dds), 1, gm.mean)
spec_dds = estimateSizeFactors(spec_dds, geoMeans = geo.mean)
spec_dds <-DESeq(spec_dds, test = "Wald", fitType = "parametric")
```

##### 6.3. 可视化 DESeq 结果：瀑布图和火山图

```r
library(ggrepel, verbose = FALSE)
#Waterfall plot for phylum

sigtabspec = results(spec_dds, contrast = c("group", "Mock", "A")) %>% 
              data.frame() %>% 
              filter(padj <=0.05, abs(log2FoldChange) > 7 ) 
sigtabspec$Order <- tax_table(spec)[,4][rownames(tax_table(spec)) %in% rownames(sigtabspec)] 
sigtabspec$Species <- rownames(sigtabspec)
#plotting for the phylum alone

spec_col=c( "deepskyblue",'springgreen3','snow3','burlywood4', 'cadetblue', 'darkblue',
          'cornflowerblue','deeppink2','orangered', 'dimgrey', 'red','limegreen',
        'cyan1','darkmagenta', 'purple', 'cyan4', 'gold', '#470e19', '#124435', '#1d723f', '#57e9ff') 


colindex = data.frame(color = spec_col[1:length(unique(sigtabspec$Order))], 
                      Order = sort(unique(sigtabspec$Order)))

ords = unique(data.frame(sigtabspec$Order)) %>% pull
colors = c()
for(i in ords){
   colors[i] = colindex[colindex$Order == i,1]
}

# Phylum order
x = tapply(sigtabspec$log2FoldChange, sigtabspec$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabspec$Order = factor(as.character(sigtabspec$Order), levels=names(x))

#Species reorder
x = tapply(sigtabspec$log2FoldChange, sigtabspec$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabspec$Species = factor(as.character(sigtabspec$Species), levels=names(x))


waterfall_p = ggplot(sigtabspec, aes(y=Species, x=log2FoldChange), stroke = 0.5) + 
  geom_vline(xintercept = 0.0, 
             color = "orange", size = 0.5, lty = 2) +
  geom_point(aes(fill = Order), alpha = 0.6, size = 10, 
  color = "black", shape = 21, stroke = 0.5) + theme_bw() + 
           scale_x_continuous(limits = c(-20, 15), n.breaks = 10) + 
           ggtitle("Log2FoldChange of Species", "Mock vs. A") +  
        scale_fill_manual(values = colors[names(colors) %in% sigtabspec$Order]) 
           
#ggsave("./deseq2/mucus/difabund_muc_DiarNoInfl_vs_NoDiar.jpeg", device = "jpeg", dpi = 300)

waterfall_p
```
```r
rm(sigtabphyl, alpha, x, colors, colindex, phyla)
```
```r
#Volcano plot

alpha = 0.05
spec_dat = results(spec_dds, contrast =  c("group", "Mock", "A"))%>% data.frame

spec_dat = spec_dat[complete.cases(spec_dat),]

spec_dat$Significant = ifelse(spec_dat$padj <= alpha, paste0("FDR < ", alpha), "Not Sig") %>% 
factor(levels = c("FDR < 0.05", "Not Sig"))
 
spec_taxa = tax_table(spec) %>% as.matrix

sigtabspec = cbind(as(spec_dat, "data.frame"), 
                   as(spec_taxa[rownames(spec_taxa) %in% rownames(spec_dat),], "matrix"))


#a costumized color scheme
phylcol=c('coral4', "deeppink",'brown2','antiquewhite4', 'cornflowerblue', 'plum4',
          'darkgoldenrod3','aquamarine4', 'yellow', 'red', 'darkblue', 'Maroon', 'Gray',
        'steelblue2','darkgreen', 'tomato1', 'cyan4', 'magenta')

colindex = data.frame(color = phylcol[1:length(unique(tax_table(spec)[,2]))], phylum = sort(unique(tax_table(spec)[,2])))


phyla = unique(data.frame(tax_table(spec)[,2])) %>% pull
colors = c()
for(i in phyla){
   colors[i] = colindex[colindex$Phylum == i,1]
}



#filtering out the taxa below 2 LFC
#sigtabspec = sigtabspec[abs(sigtabspec$log2FoldChange)>2,]

volc_p = sigtabspec %>% group_by(log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% 
ggplot(aes(x = log2FoldChange, y = -log10(pvalue), label = Genus)) + 
geom_hline(yintercept = -log10(sigtabspec[sigtabspec$Significant == "Not Sig","pvalue"])  %>% 
           max, color = alpha("red",0.5), lty = 2) +
geom_vline(xintercept = 0, color = alpha("black", 0.3))  +
geom_point(data = sigtabspec[sigtabspec$Significant == "Not Sig",], 
           aes(x = log2FoldChange, y = -log10(pvalue)), 
           color = alpha("darkgreen", 0.6),  size = 2)  + 
theme_bw(base_size = 12) + 
theme(legend.position= "right", 
      text = element_text(size = 15, face = "bold")) + 
geom_point(data = sigtabspec[sigtabspec$Significant == "FDR < 0.05",],
    aes(x = log2FoldChange, y = -log10(pvalue), 
        fill = Phylum), size = 6, alpha = 0.5,  
        color = "black", shape = 21, stroke = 0.5) +  
scale_fill_manual(values = colors[names(colors) %in% 
                sigtabspec[sigtabspec$Significant == "FDR < 0.05", "Phylum"]]) +
geom_text_repel( nudge_y = 0.15, nudge_x = -.5, 
               data= top_n(sigtabspec[sigtabspec$Significant == "FDR < 0.05" & sigtabspec$log2FoldChange < -2,], 
               n =  -10, wt = pvalue),
               aes(label = Genus), 
               size = 2.5, 
               box.padding = unit( 0.4, units ="lines"),
               point.padding = unit(0.4, "lines"), max.overlaps = 20)  + 
geom_text_repel(nudge_y = 0, nudge_x =0.5, 
               data= top_n(sigtabspec[sigtabspec$Significant == "FDR < 0.05" & sigtabspec$log2FoldChange > 2,],
               -10, pvalue), 
               aes(label = Genus), size = 2, 
               box.padding = unit( 0.4, units ="lines"),
               point.padding = unit(0.4, "lines"), max.overlaps = 20)+
geom_text( aes(x = 4, y =0, label = "Not Sig"), color = "red", size = 2.5) +
geom_text(aes(x = 4, y =5, label = "FDR < 0.05"), color = "red", size = 2.5)+ 
ggtitle (label = "Volcano Plot of the top 10 most significant log2FoldChange", "Species in Mock vs. A groups") + 
scale_y_continuous(limits = c(0, 10), n.breaks = 5) + 
scale_x_continuous(limits = c(-7.5, 6), n.breaks = 10) + 
guides(size = "none") 

volc_p
```
```r
#ggsave("./deseq2/mucus/volc_gen_DiarNoInfl_vs_NoDiar_muc.jpeg", device = "jpeg", dpi = 300)
```

#### 7. 生物标志物和分类单元数据之间的热图关联

```r
#correlation heatmap between biomarker and the Genus
gen <- gloomer(ps_rar, taxa_level = "Genus", NArm = TRUE)

gen_log = filter_taxa(gen, function(x) sum(x>0)>0, TRUE) #filtering the zero counts out
gen_log = transform_sample_counts(gen_log, function(x) {log(1+x)})

# Creating dummy variables                                       
biodf = sample_data(gen_log)

biodf = data.frame(biodf, geneA = rnorm(87, mean = 0, sd = 1), 
          geneB = rnorm(87, mean = 0, sd = 0.5), 
          geneC = rnorm(87, mean = 0, sd = 0.005)) %>% 
          select(-c(1:2)) %>% 
          as.matrix()

                     
countdf = as.matrix(otu_table(gen_log))       
 
countdf = countdf[!rownames(countdf) %in% c("Unknown", "uncultured", "Uncultured", "Unassigned", "NA"),] 
                                  
countdf = t(countdf)  
countdf = countdf[rownames(countdf) %in% rownames(biodf),] 
        

#cor gene taxa

cor_main = Hmisc::rcorr(countdf, biodf, type = "spearman")
cors = cor_main$r
cors = cors[rownames(cors) %in% colnames(countdf), colnames(cors) %in% colnames(biodf)] #rows as taxa, cols as chemical data

#calculating the qvalues for the correlations
cor.pval = cor_main$P[rownames(cor_main$P) %in% rownames(cors) , colnames(cor_main$P) %in% colnames(cors)] 
cor.qval = p.adjust(cor.pval, method = "BH")

q.vals = matrix(cor.qval, ncol = ncol(cor.pval), nrow = nrow(cor.pval), dimnames = list(rownames(cor.pval), colnames(cor.pval)))
                     
#Adding significance signes to the qval matrix to be used in the heatmap later on
q.vals[cor.qval <0.05] = "*"
q.vals[cor.qval >= 0.05] = ""

library(Biobase, verbose = FALSE)
asvdat = as(cors,"matrix")#in the aassay dataset, we add our correlation matrix instead of the abundance matrix

taxadat = Biobase::AnnotatedDataFrame(data.frame(cors))#taxa table
pdata = cors %>% data.frame

x = ExpressionSet(assayData = asvdat,  featureData = taxadat  )

#Adding phenotype data
pData(x) <- pdata          

# Filtering based on row standard deviation and choosing the most variable 50 taxa
library(matrixStats)
sds <- rowSds(Biobase::exprs(x))
o <- order(sds, decreasing = TRUE)[1:50]
h_1 <- hclust(dist(Biobase::exprs(x)[o,]), method = "ward.D2")
h_2 <- hclust(dist(t(Biobase::exprs(x)[o,])), method = "ward.D2")

#making a phylum annotation and it only accepts one column dataframe

row.annot = gen_log@tax_table[rownames(gen_log@tax_table) %in% rownames(Biobase::exprs(x)[o,]),2]

#making color index for the phylum annotation 

phylcol=c('coral4', 'cyan','#ff00aa', 'tomato1', 'cornflowerblue', 'plum4',
          'darkgoldenrod3','aquamarine4', 'cadetblue2', 'red', 'darkblue', 'Maroon', 'Gray',
        'steelblue2','darkmagenta', 'antiquewhite4', "darkorange", 'darkgreen')
set.seed(2)
phylcol = sample(phylcol, size = length(unique(row.annot[,1])), replace = F)
phyl.col = data.frame(Phylum = unique(row.annot[,1]), phyl.col = phylcol[1:length(unique(row.annot[,1]))])
rownames(phyl.col) <- NULL
phyl.col = column_to_rownames(phyl.col, "Phylum") %>% as.matrix

library(RColorBrewer, verbose = FALSE)


 
mat = matrix(NA, ncol = 1, nrow = nrow(row.annot))
for(i in 1:nrow(row.annot)){

 mat[i,] = ifelse(row.annot[i,1][[1]] %in% rownames(phyl.col), phyl.col[rownames(phyl.col) %in% row.annot[i,1][[1]] ,1], "NA")

}
colnames(mat) <- "col"

row.annot = cbind(row.annot, mat) %>% data.frame()

pheatmap(Biobase::exprs(x)[o,], annotation_row = row.annot %>% select(1),
                          cellheight = 12, annotation_colors = list(
                            Phylum = phyl.col[,1]),
                           cellwidth = 15, cutree_rows = 4, border_color = NA, 
                           fontsize_number = 15, number_color = "black", 
                           display_numbers = q.vals[rownames(q.vals) %in% 
                           rownames(Biobase::exprs(x)[o,]), colnames(q.vals) %in% 
                           colnames(exprs(x)[o,])], angle_col = 45,
                           Rowv = as.dendrogram(h_1), Colv = as.dendrogram(h_2), 
                           cutcluster_rows = T, cluster_cols = F, col = brewer.pal(9, "Reds"), 
                           width = 10, height = 12, main = "Spearman correlation of top 50 Genre \nand 3 Genes, clustered row ")
```
```r
#ggsave(plot = pheat.chem.scfa, "./heatmap/heatmap.SCFA.gen_50.jpeg", dpi = 750, height = 10, width = 8)
```
在原假设下，p 值是具有均匀分布的随机值。因此，我们必须执行调整后的 p 值来校正族错误率 (FWER)。FWER 是指同时进行多个假设检验时至少出现一个 I 类错误（误报）的概率。

```r
# Set the parameters
num_simulations <- 10000  # Number of simulated datasets
sample_size <- 50         # Sample size for each group

# Initialize a vector to store p-values
simulated_p_values <- numeric(num_simulations)

# Simulate data and compute p-values
for (i in 1:num_simulations) {
  # Simulate data under the null hypothesis for two groups
  group1 <- rnorm(sample_size)
  group2 <- rnorm(sample_size)
  
  # Perform a two-sample t-test on the simulated data
  test_result <- t.test(group1, group2)
  
  # Store the p-value
  simulated_p_values[i] <- test_result$p.value
}

# Plot a histogram of simulated p-values
hist(simulated_p_values, breaks = 30, col = "lightblue", main = "Distribution of Simulated P-values under Null hypothesis")
```
```r
print(glue("{sum(simulated_p_values<=0.05)} of the tests out of 10000 test were significant. \n This indicates that under the Null hypothesis alpha % of times your test will be rendered false positive (type-I error)."))
```
```r
## 515 of the tests out of 10000 test were significant. 
## This indicates that under the Null hypothesis alpha % of times your test will be rendered false positive (type-I error).
```

#### 8. 系统发育树可视化

```r
library(phytools)
library(TDbook)
library(ggimage)
library(treeio)
library(tidyverse)
library(tidytree)


spec <- gloomer(ps_rar, taxa_level = "Class")
spec <- transform_sample_counts(spec, function(x){x/sum(x)*100})
spec <- prune_taxa(taxa_sums(spec)>1, spec)

taxadf <- as.data.frame(tax_table(spec))

color.index <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")



#creating a highlihgt df with a costum function


mrca.wrap <- function(highlight, taxa.df, tree, tax_level="Phylum", type.id = "node"){

taxa <- taxa.df %>% data.frame 
node.id <- list()
    
for(i in phyls){
  node.id[i] <- findMRCA(tree, tree$tip.label[taxa[,tax_level] == i], type = type.id) 
    node.df <- data.frame(phyl = names(node.id), node.id = c(node.id[1][[1]]))
}
    node.id = node.id = as(node.id,"matrix")
    node.df = data.frame(highlights = rownames(node.id), node.id = unlist(node.id))
    return(node.df)
    }

ps_mock = phyloseq::subset_samples(spec, group == "Mock") 
phyls = tax_table(ps_mock)[rowSums(ps_mock@otu_table)>0,2] %>% data.frame() %>% distinct() %>% pull()#chosoing phylum that that are in Mock samples

node.df <- mrca.wrap(highlight = phyls, 
          taxa.df = tax_table(spec), 
          tree = phy_tree(spec), tax_level = "Phylum", 
          type = "node")
          

#drwing the tree
tree_p =  ggtree(spec, aes(color = Phylum),
branch.length = "none",
layout = "circular",
show.legend =TRUE,
open.angle = 5,
size = 1.5)  +  
ggtitle("Phylogenetic tree of different Classes.", 
      "Mock comunity phyla are highligted!") +
geom_tiplab(aes(label=Class, color = Phylum), 
check.overlap = FALSE, 
face = "bold", 
size = 3, 
offset = 0.3
)  +
scale_color_manual(values = sample(color.index, 
                    length(unique(taxadf$Phylum)), FALSE), 
                    breaks = unique(taxadf$Phylum)) + #to remove the NA from the legend
geom_highlight(data = node.df, 
lwd = 0.25, lty = 3,  alpha = 0.1,
aes(node = node.id, fill = highlights), 
extend =0.05, to.bottom = TRUE, 
align = "right", show.legend = TRUE) +
geom_nodepoint( pch = 21, color = "black", size = 2, fill = "white") + 
geom_nodepoint(aes(subset = node %in% node.df$node.id), pch = 21, size = 3, color = "black", 
fill = c("#00ffff", "#ff7700")
) +
geom_tippoint(size = 1, pch = 21,  color = alpha(colour = "white", alpha = 0.2)) + 
geom_label(aes(x = branch, label = round(branch.length,2)),
label.padding = unit(0.05, "line"), size = 3, inherit.aes = TRUE) +
scale_fill_manual(values = c( "#16f476", "#fc00fc")) 


# geom_cladelab(node = node.df$node.id, align = FALSE,
# label = node.df$highlights,
#  offset.text = 0.25, 
#  offset = 16,
#  barsize = 0.05, 
#  extend = 0.01,
#  angle ="auto",
#  fontface = 2, size = 7) +



#Plotting abundance on tree
rel.df <- psmelt(spec) %>% select(OTU, Abundance, group) %>% group_by( OTU, group) %>% summarise(rel = mean(Abundance), .groups = "drop")  %>% pivot_wider(values_from = "rel", id_cols = "OTU", names_from = "group") %>% column_to_rownames( "OTU")
rescaled.reldf <- apply(rel.df, 2, function(x) {log(x+1)})


library(ggnewscale, verbose = FALSE)
p2 <- tree_p + new_scale_fill()
p3 <- gheatmap(p2, font.size = 5,  rescaled.reldf, offset=6.04, width=.8, 
         colnames_angle=60, colnames_offset_y = -.45, colnames_position = "top") +
    scale_fill_viridis_c(option="C", name="Log relative abundance, %")

p3
```

原文链接和[Github存储库位置](https://github.com/farhadm1990/MAC2023.github.io)