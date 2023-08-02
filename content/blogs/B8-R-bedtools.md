---
title: 使用R语言实现bedtools求交集的功能？
date: 15th Feb 2023
description: 使用R语言实现bedtools求交集的功能.
image: /blogs-img/R.webp
alt: 使用R语言实现bedtools求交集的功能？
ogImage: /blog-img/R.webp
tags: ['r语言']
published: true
---

Bedtools作为基因组研究的 “ 瑞士军刀 ”， 功能强大且易于操作，是生信行业不可多得的好软件。通常对bed区间的注释，我们使用其中“ 求交集 ”的功能（bedtools intersect) ，但是有一个很不方便的地方，我们通常要生成对应的bed文件，再注释完成后还需要用R语言等读入才能继续分析，所以整合度不是很好，本文希望提供R语言的思路来解决该问题。


### 什么是bedtools?
### **bedtools**：*一个强大的基因组算法工具集*

总的来说，**bedtools**实用程序是用于广泛基因组学分析任务的瑞士军刀。最广泛使用的工具支持*基因组算法*：即基因组的集合论。例如，**bedtools**允许从广泛使用的基因组文件格式（如 BAM、BED、GFF/GTF、VCF）的多个文件中*交叉*、*合并*、*计数*、*补充*和*洗牌基因组区间。*虽然每个单独的工具都旨在完成一项相对简单的任务（例如， *相交*两个间隔文件），可以通过在 UNIX 命令行上组合多个 bedtools 操作来进行相当复杂的分析。

**bedtools**是在犹他大学的昆兰实验室开发的，受益于世界各地科学家的杰出贡献。你可以从该链接[https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)获取它所有的所有功能。

### 第一个问题：获取bed区间开始/结束位点落在基因区间上的结果

这个问题用bedtools也不好解决，因为实际上我们需要获取交集的子集。需要首先把bed区间转变为只有1bp的区间，这增加了操作步骤，但是在R语言中我们可以很灵活实现该方案。

需要先导入R包


```r
library(GenomicRanges)
library(readr)
library(dplyr)
```

然后写一个函数，负责将bed文件转化为GRanges对象。

```R
makeGranges<-function(x,meta.name){
  df<-read_tsv(x,col_names=F)
  if(stringr::str_detect(df$X1[1],"chr")){
    df$X1<-stringr::str_remove(df$X1,"chr")
  }
  names(df)<-c("seqname","start","end",meta.name)
  makeGRangesFromDataFrame(df,keep.extra.columns = T,ignore.strand = T)
}
```
然后就是实现该题目功能的函数了

```r
anno_start<-function(x,y){
  x_tmp<-narrow(x,start=1,width = 1)
  res<-findOverlaps(x_tmp,y,ignore.strand=T)
  jiaoji<-x[queryHits(res)]
  mcols(jiaoji)$anno<-mcols(y[subjectHits(res)])[[1]] ## assume that the first meta column should be anno info
  jiaoji
}
```
其实就是通过缩短bed区间到起始位置1bp，然后求交集即可。

### 第二个问题：如何实现百分比交集

我们通常不仅仅想知道交集，还想知道交集之间相交部分 的长度占自身的百分比来控制结果的输出。

```r
get_overlap_percentage<-function(x,y,query.restrict=T,pct=0.2){
  hits <- findOverlaps(x,y,ignore.strand=T)
  ints<-pintersect(x[queryHits(hits)],y[subjectHits(hits)])
  if(query.restrict==TRUE){
    int_p<-width(ints)/width(x[queryHits(hits)])
  }else{
    int_p<-width(ints)/width(y[subjectHits(hits)])
  }
  hits <- hits[int_p>=pct]
  jiaoji<-x[queryHits(hits)]
  mcols(jiaoji)$anno<-mcols(y[subjectHits(hits)])[[1]]
  jiaoji
}
```
虽然也不难，但是这个功能实现起来，还是比第一个问题要复杂的，因为交集会存在一对多的情况，所以我们需要两次求交，分别计算长度以及百分比，然后控制第一次的输出。