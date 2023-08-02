---
title: 基于Julia语言的多线程barcode 拆分
date: 8th May 2021
description: 基于Julia语言的多线程barcode 拆分.
image: /blogs-img/julia.webp
alt: 基于Julia语言的多线程barcode 拆分
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

### 拆分原理

- 软件的逻辑是首先获取barcode列表。然后采用多线程分别在fastq文件中并行提取对应barcode的reads。
 
- WGS的下机数据经常出现在fastq2里。所以程序会从fastq中自动查找是否存在对应barcode。
 
- 程序可以自动检测barcode始于开始还是末尾，计算hanming距离，运行1bp的mismatch。 
 
### 方法


```bash
julia -t number threads split_barcode.jl -b barcode list -r read2 file -l read1 file -o outpath
```
 
### 代码举例



```julia
using BioSequences
using FASTX
using CodecZlib
using StringDistances
using Base.Threads
using Getopt

function read_barcode(barcode_fs::String)
    labels = String[]
    barcodes = String[]
    if isfile(barcode_fs)
        for line in eachline(barcode_fs)
            label, sequence = split(line)
            push!(labels, label)
            push!(barcodes, sequence)
        end
    else
        println("no barcode file detected!")
    end
    labels, barcodes
end

function detective_barcode(seq::LongDNASeq, barcode::LongDNASeq)
    barcode_l = length(barcode)
    starter = seq[1:barcode_l]
    ender = seq[end-barcode_l:end]
    dist = Hamming()(starter, barcode)
    if dist == 0 || dist == 1 ## if barcode in title
        return 1, 1
    else
        dist = Hamming()(ender, barcode)# calular hamming distance allow 1bp mismatch
        if dist == 0 || dist == 1
            return 1, 2
        else
            return -1, 0
        end
    end
end

function process(fq_file_1::String, fq_file_2::String, barcode::LongDNASeq, label::String, outpath::String)
    reader_2 = FASTQ.Reader(GzipDecompressorStream(open(fq_file_2)))
    reader_1 = FASTQ.Reader(GzipDecompressorStream(open(fq_file_1)))
    wt_1 = FASTQ.Writer(open("$(outpath)/barcode_$(label)_1.fq", lock=true, "a")) ## lock for multiple threads safety
    wt_2 = FASTQ.Writer(open("$(outpath)/barcode_$(label)_2.fq", lock=true, "a"))
    for (record_1, record_2) in zip(reader_1, reader_2)
        fm, weizhi = detective_barcode(FASTQ.sequence(record_2), barcode)
        if fm == 1
            write(wt_1, record_1)
            if weizhi == 1
                write(wt_2, FASTQ.Record(FASTQ.identifier(record_2), FASTQ.sequence(record_2)[11:end], FASTQ.quality(record_2)[11:end]))
            elseif weizhi == 2
                write(wt_2, FASTQ.Record(FASTQ.identifier(record_2), FASTQ.sequence(record_2)[1:100], FASTQ.quality(record_2)[1:100]))
            end
        end
    end
    close(reader_1)
    close(reader_2)
    close(wt_1)
    close(wt_2)
end


function split_barcode(barcodes_fs::String, fq_file_2::String, fq_file_1::String, outpath::String)
    labels, barcodes = read_barcode(barcodes_fs)
    barcodes = LongDNASeq.(barcodes)
    duiying = Dict(zip(barcodes, labels))
    @threads for barcode in barcodes
	    label = duiying[barcode]
        process(fq_file_1, fq_file_2, barcode, label, outpath)
    end
end

function Argparse()
    lst = Dict{String,String}()
    for (opt, arg) in getopt(ARGS, "hb:r:l:o:", ["help","barcodes=", "read2=","read1=","output="])
        opt = replace(opt, "-" => "")
        arg = replace(arg, "=" => "")
        if opt == "help" || opt == "h"
            println("this program is for spliting fastq file into different part according barcodes!")
            println("-h\t--help\tshow this message")
            println("-b\t--barcodes\tinput your barcodes file, should be negetive sequennce!")
            println("-r2\t--read2\tinput your fastq file read_2")
            println("-r1\t--read1\tinput your fastq file read_1")
            println("-o\t--output\tinput output file path")
        elseif opt == "barcodes" || opt == "b"
            lst["barcodes"] = arg
        elseif opt == "read2" || opt == "r"
            lst["read2"] = arg
        elseif opt == "read1" || opt == "l"
            lst["read1"] = arg
        elseif opt == "output" || opt == "o"
            lst["output"] = arg
        else
            println("please check your parameter!")
        end
    end
    lst
end


args = Argparse()
if length(args) > 1
    split_barcode(args["barcodes"], args["read2"], args["read1"], args["output"])
end
```