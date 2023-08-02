---
title: Julia语言模仿BAM文件的pileup类似操作
date: 22th Mar 2023
description: Julia语言模仿BAM文件的pileup类似操作.
image: /blogs-img/julia.webp
alt: Julia语言模仿BAM文件的pileup类似操作
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

### 简介

使用过`pysam`和`samtools`的小伙伴肯定了解 `pileup`的操作，如果把BAM文件看作表格的话，那么通常我们是按行去解析它的record，进而获得一些信息，例如比对到哪条染色体，比对的开始位置和结束位置等. 另一种情况下，我们想要按照列去循环解析，得到这个列上的具体信息，典型的就是这个列上比对序列的碱基是什么？比对序列的位置是什么？以及是Match or Mismatch or indel 等。那么，该操作就需要引入`pileup`操作了。

`pysam`包中分别有`pileup`，`PileupColumn`以及`PuleupRead` 对象来承担上述任务, 那么我们简单的在julia使用，不需要构建对象这么复杂的操作，写几个函数即可
### 引入需要的包

```r
using XAM
using BioAlignments
using GenomicFeatures
using BioGenerics
```

### 获取比对到某个位点和某个区间的read 数量

```r
function pileup(bam::BAM.Reader,contig::String, pos::Int64)
    ### get the intervealCollection
    site = Interval(contig, pos,pos)
    ### get the pileup
    GenomicFeatures.eachoverlap(bam,site)
end

## return number of segments in the pileup
function nsegments(bam::BAM.Reader,contig::String, start::Int64, stop::Int64)
    interval = Interval(contig, start, stop)
    num = 0
    for i  in eachoverlap(bam,interval)
        num+=1
    end
    num
end

function nsegments(bam::BAM.Reader,contig::String, pos::Int64)
    interval = Interval(contig, pos, pos)
    num = 0
    for i  in eachoverlap(bam,interval)
        num+=1
    end
    num
end
```
### 获取比对到某个位点的Query 序列的位置和比对情况

```r
function get_query_position(bam::BAM.Reader, contig::String, pos::Int)
    # get the pileup at the contig position
    PileupReads = pileup(bam, contig, pos)
    # get the reference base
    Iterators.map(x->(ref2seq(BAM.alignment(x),pos)[1]),PileupReads)
end
    
function get_query_operation(bam::BAM.Reader,contig::String, pos::Int)
    # Get pileup reads
    PileupReads = pileup(bam,contig, pos)
    # Map each read to the base at the given position
    Iterators.map(x->(ref2seq(BAM.alignment(x),pos)[2]),PileupReads)
end
```

该`get_query_position`函数返回一个生成器，可以迭代获取比对到某个参考位点的record上的具体位置，例如

```r
bam = open(BAM.Reader,"45T.cs.rmdup.sort.bam",index="45T.cs.rmdup.sort.bam.bai")
pos = get_query_position(bam,"chr6",36434531)
for p in pos
    println(p)
end
```
当然，你可以做一个双重循环，计算很多位点的具体位置。

同样，该`get_query_operation`函数返回的是具体比对类型，我们都知道，比对有M(match),I(insertion),D(deletion)组成，这个就可以打印对应的行为。

```r
bam = open(BAM.Reader,"45T.cs.rmdup.sort.bam",index="45T.cs.rmdup.sort.bam.bai")
pos = get_query_operation(bam,"chr6",36434531)
for p in pos
    println(p)
end
M
M
M
M
M
...
```