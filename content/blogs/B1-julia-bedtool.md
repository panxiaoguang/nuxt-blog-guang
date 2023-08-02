---
title: 用julia实现bedtools的intersect-bed功能
date: 7th Oct 2020
description: 用julia实现bedtools的intersect-bed功能.
image: /blogs-img/julia.webp
alt: 用julia实现bedtools的intersect-bed功能
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

### 前言
自己写的好几种算法企图实现bedtools的功能，虽然julia性能足够好，但都难以在效率上达到bedtools的性能，于是最后只能借助轮子了。

### 代码

```julia
using CSV
using DataFrames
using Tables
using GenomicFeatures

function fs2NameTuple(fs::String)
    query = DataFrame(CSV.File(fs, delim="\t", header=0))
    if size(query)[2] == 5
        query = rename!(query, :Column1 => :Chromosome, :Column2 => :Start, :Column3 => :End, :Column4 => :Name, :Column5 => :Score)
        # query = sort!(query, [:Chromosome,:Start])
        query = Tables.rowtable(query)
        result = IntervalCollection([Interval(qy.Chromosome, qy.Start, qy.End, '?', qy.Name) for qy in query], true)
        result
    elseif size(query)[2] == 4
        query = rename!(query, :Column1 => :Chromosome, :Column2 => :Start, :Column3 => :End, :Column4 => :Name)
        # query = sort!(query, [:Chromosome,:Start])
        query = Tables.rowtable(query)
        result = IntervalCollection([Interval(qy.Chromosome, qy.Start, qy.End, '?', qy.Name) for qy in query], true)
        result
    else
        println("please confirm your bed files")
    end
end

function getInterval(A::IntervalCollection, B::IntervalCollection)
    for i in eachoverlap(A, B)
        println(i[1].seqname, "\t", i[1].first, "\t", i[1].last, "\t", i[1].metadata, "\t", i[2].seqname, "\t", i[2].first, "\t", i[2].last, "\t", i[2].metadata)
    end
end
getInterval(fs2NameTuple("test1.bed"),fs2NameTuple("test2.bed"))
```

### 原理

首先通过`CSV`将bed文件读入内存，然后将`dataFrame`转化为`nametuple`，最后在转化为`intervalCollection`。
通过内置的`overlap`功能可以几乎瞬间获取所有的交集，性能已经达到了`c`语言的程度。