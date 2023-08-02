---
title: 一文读懂lncRNA分析
date: 17th Feb 2020
description: 一文读懂lncRNA分析.
image: /blogs-img/lncRNA.webp
alt: 一文读懂lncRNA分析
ogImage: /blogs-img/lncRNA.webp
tags: ['转录组', '生物信息', 'r语言']
published: true
---

### 1:对测序下机数据进行质量检测（QC）

使用软件：**fastqc**

命令行参数：

```bash
-o --outdir:输出路径
--extract：结果文件解压缩
--noextract：结果文件压缩
-f --format:输入文件格式.支持bam,sam,fastq文件格式
-t --threads:线程数
-c --contaminants：制定污染序列。文件格式 name[tab]sequence
-a --adapters：指定接头序列。文件格式name[tab]sequence
-k --kmers：指定kmers长度（2-10bp,默认7bp）
-q --quiet： 安静模式
```


运行命令：

```bash
fastqc -f fastq -t 8 -o ./QCreport/ fastq1 fastq2
```

注意：当样本数量很多时，我们可以采用 multiqc 软件将 fastqc 软件的结果进行合并

```shell
multiqc ./*.zip
```


### 2.对 qc 后的数据进行过滤

使用软件：**fastp**

命令行参数：

```bash
usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
options:
  # I/O options   即输入输出文件设置
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                         read2 output file name (string [=])
  -6, --phred64                      indicates the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2. (int [=2])
    --reads_to_process               specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])

  # adapter trimming options   过滤序列接头参数设置
  -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])

  # global trimming options   剪除序列起始和末端的低质量碱基数量参数
  -f, --trim_front1                  trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                   trimming how many bases in tail for read1, default is 0 (int [=0])
  -F, --trim_front2                  trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                   trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])

  # polyG tail trimming, useful for NextSeq/NovaSeq data   polyG剪裁
  -g, --trim_poly_g                  force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
      --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
  -G, --disable_trim_poly_g          disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data

  # polyX tail trimming
  -x, --trim_poly_x                    enable polyX trimming in 3' ends.
      --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])

  # per read cutting by quality options   划窗裁剪
  -5, --cut_by_quality5              enable per read cutting by quality in front (5'), default is disabled (WARNING: this will interfere deduplication for both PE/SE data)
  -3, --cut_by_quality3              enable per read cutting by quality in tail (3'), default is disabled (WARNING: this will interfere deduplication for SE data)
  -W, --cut_window_size              the size of the sliding window for sliding window trimming, default is 4 (int [=4])
  -M, --cut_mean_quality             the bases in the sliding window with mean quality below cutting_quality will be cut, default is Q20 (int [=20])

  # quality filtering options   根据碱基质量来过滤序列
  -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])

  # length filtering options   根据序列长度来过滤序列
  -L, --disable_length_filtering     length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])

  # low complexity filtering
  -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
  -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])

  # filter reads with unwanted indexes (to remove possible contamination)
      --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

  # base correction by overlap analysis options   通过overlap来校正碱基
  -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled

  # UMI processing
  -U, --umi                          enable unique molecular identifer (UMI) preprocessing
      --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])

  # overrepresented sequence analysis
  -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
  -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])

  # reporting options
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])

  # threading options   设置线程数
  -w, --thread                       worker thread number, default is 3 (int [=3])

  # output splitting options
  -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
  -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
  -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])

  # help
  -?, --help                         print this message
```

运行命令：

```bash
fastp -c -i fastq1 -o ./cleandata/ -I fastq2 -O ./cleandata/
```

### 3.对干净数据去除 rRNA 操作（存疑）

在 ncbi 下载 rRNA 数据保存为 hg19_rRNA.fasta

使用 bwa 建立索引文件

```bash
bwa index hg19_rRNA.fasta
```

使用 bwa 将测序数据比对到 rRNA 数据库，去除比对上的数据，并恢复为 fastq 文件

```bash
bwa mem hg19_rRNA.fasta fastq1 fastq2 | samtools view -bf4 |samtools -c 9 -1 paired1.fq.gz -2 paired2.fq.gz -
```

### 4.将去除 rRNA 的数据比对到 hg19|hg38 上

**使用软件**：subread

**建立索引**：

```bash
subread-buildindex -o hg19 ./database/hg19.fasta
```

**比对**：

```bash
subjunc -T 5 -i hg19 -r fastq1 -R fastq2 -o out.bam
```

**统计 readcounts**

```bash
featureCounts -p -T 6 -t exon -g transcript_id -a all_lncRNA_know.gtf -o sample.lncRNA.readcount out.bam
```

> 这里需要注意：因为一个基因会包含多个转录本，一个转录本又包含多个外显子和内含子，所以为了精确比对结果，我们需要基于外显子进行比对`-t exon`，软件会自动根据注释文件进行合并同一个转录本的计数结果，并计算出转录本的长度汇总到结果中.

> 另外，需要注意的是，正因为可变剪切的存在，我们不能简单的认为转录本的计数结果相加即为对应基因的计数结果，这是不准确的，所以，要想知道基因的计数结果，我们还需要修改参数，另行计算`-g gene_id`

**将计数的结果提取出来**：

```bash
awk -F '\t' '{print $1,$7}' OFS='\t' sample.lncRNA.readcount > ./tiqu/sample.counts
```

> 注：对于结果中的 length，有人以为就是有效长度（effctive length),但事实是，这个长度仅仅只是一个转录本上所有外显子长度的加和而已，因此，并不能作为有效长度。网上有估算有效长度的公式：

> ```r
> efflength = transcripts_length - LFD +1
> ```
>
> LFD 是测序时的平均片段长度，这个可以估算，但不准确，最好的方法是询问测序人员这个数值。

### 5.整合所有样本信息到表达矩阵并计算 TPM

计算完 readcount，转录组测序的打基础过程已经完毕，后续的分析都是基于此的，但是，为了下游分析的方便，我们有必要将所有的数据整理到一个表达矩阵中。

使用**python**或者**R**或者**shell**，我建议使用 python,因为简单快速易控：

```python
import pandas as pd
import glob
#使用pandas速度要更快，先预读一个样本作为dataframe的起始框架
df=pd.read_csv("sample.lncRNA.trans.counts",sep="\t",skiprows=1)
#featurecount的计算结果提取后通常有这样的表头
df=df[['Geneid','LncRNA/bidui/normal/sample.subjunc.bam']]
df.columns=['id','counts']
#遍历所有样本，使用pd.join整合到一个大表
for name in glob.glob("*.counts"):
    sample=name.replace(".lncRNA.trans.counts","")
    df2=pd.read_csv(name,sep="\t",skiprows=1)
    df2.columns=['id','counts','ot']
    df2=df2[['counts']]
    df2.columns=[sample]
    df=pd.concat([df,df2],axis=1)
#删除预读样本信息，并写入文件
df.drop('counts',axis=1).set_index("id").to_csv("normal_lncRNA_trans_allcounts.csv",sep="\t")
```

目前主流的表达量标准化方法多种多样，有 RPKM，FPKM，TPM，CPM 等，下游差异分析软件 deseq2 也有自己的标准化方法，但是，为了保证生物学意义，同时又保证数据的可比性，我还是建议选择 TPM。

##### reads Count

**定义**: 高通量测序中比对到 exon 上的 reads 数。可使用 featureCount 等软件进行计算。

**优点：**可有效说明该区域是否真的有表达及真实的表达丰度。能够近似呈现真实的表达情况。有利于实验验证。

**缺点：**由于 exon 长度不同，难以进行不同 exon 丰度比较；由于测序总数不同，难以对不同测序样本间进行比较。

##### RPKM/FPKM

**定义：**
RPKM: Reads Per Kilobase of exon model per Million mapped reads (每千个碱基的转录每百万映射读取的 reads)；
FPKM: Fragments Per Kilobase of exon model per Million mapped fragments(每千个碱基的转录每百万映射读取的 fragments)
**公式：**

$$RPKM =\frac{ExonMappedReads  * 10^9 }{TotalMappedReads * ExonLength}$$

上述公式可从下面公式推导而出：
$$RPKM=\frac{ExonMappedReads/ExonLength*10^9}{TotalMappedreads/GenomeLength}$$

**解释**：ExonMappedReads 即为比对到该 exon 上的 reads count； TotalMappedReads 即为比对到基因组上所有 reads count 的总和；ExonLength 为该 Exon 的长度；GenomeLength 即为基因组全长，因为是相同基因组，所以该数值可消除。

**优点：**tophat-cufflinks 流程固定，应用范围广。理论上，可弥补 reads Count 的缺点，消除样本间和基因间差异。

**讨论：**有人说 RPKM/FPKM 标准化特别不合理，看着是个大牛 YellowTree。公式 2 中，TotalMappedReads/GenomeLength 为测序深度，ExonMappedReads / ExonLength 可以简单的认为是该 Exon 上的“测序深度”。两者相除，就得出该 Exon 依据测序深度而进行的标准化，那么因 Exon 长短、测序深度造成的样本间造成的偏差，都可以消除。因一般是相同物种，基因组一般相同，所以公式 2 换算并消去 GenomeLength，就成为公式 1 的形式了。不知道哪里错了，斗胆提出质疑：RPKM/FPKM 怎么就不能消除两种类型的 bias？不过有论文陈述说 RPKM 的结果难以消除组间测序造成的差异，可能未采用比对到基因组上所有的 reads 数，而是采用了比对到所有 Exon 的 reads 数作为 TotalMappedReads 吧。不是很确定。

**FPKM：**与 RPKM 计算过程类似。只有一点差异：RPKM 计算的是 reads，FPKM 计算的是 fragments。single-end/paired-end 测序数据均可计算 reads count，fragments count 只能通过 paired-end 测序数据计算。paired-end 测序数据时，两端的 reads 比对到相同区域，且方向相反，即计数 1 个 fragments；如果只有单端 reads 比对到该区域，则一个 reads 即计数 1 个 fragments。所以 fragments count 接近且小于 2 \* reads count。

##### RPM

**定义：**RPM/CPM: Reads/Counts of exon model per Million mapped reads (每百万映射读取的 reads)

**公式：**

$$RPM = \frac{ExonMappedReads * 10^6}{TotalMappedReads}$$

**优点：**利于进行样本间比较。根据比对到基因组上的总 reads count，进行标准化。即：不论比对到基因组上的总 reads count 是多少，都将总 reads count 标准化为 10^6。

**缺点：**未消除 exon 长度造成的表达差异，难以进行样本内 exon 差异表达的比较。

##### TPM

**定义：**TPM: Transcripts Per Kilobase of exon model per Million mapped reads (每千个碱基的转录每百万映射读取的 Transcripts)

**公式：**

$$TPM=\frac{N_i/L_i*10^6}{sum(N_1/L_1+N_2/L_2+··+N_n/L_n)}$$

**解释：**Ni 为比对到第 i 个 exon 的 reads 数； Li 为第 i 个 exon 的长度；sum(N1/L1+N2/L2 + ... + Nn/Ln)为所有 (n 个)exon 按长度进行标准化之后数值的和。

**计算过程：**首先对每个 exon 计算 Pi=Ni/Li，即按长度对 reads count 进行标准化；随后计算过程类似 RPM (将 Pi 作为正常的 ExonMappedReads，然后以 RPM 的公式计算 TPM)。

**优点：**首先消除 exon 长度造成的差异，随后消除样本间测序总 reads count 不同造成的差异。

**缺点：**因为不是采用比对到基因组上的总 reads count，所以特殊情况下不够准确。例如：某突变体对 exon 造成整体影响时，难以找出差异。

##### 相互关系

**评价：**以上几种计算 exon 表达丰度的方法，差异不是非常大。如果结果是显著的，那么采用上面任一计算方法大多均可找出显著结果。但是当表达风度差异不是那么显著时，不易区分不同类别，需要根据实际需要选择对应的标准化方法。

**注意：**以上 TotalMappedReads 推荐首选比对到基因组上的总 reads 数，而不是比对到 exon 或者 gene 上总 reads 数。这同样需要根据实际情况而确定。

##### 如何根据 counts 计算 TPM？

这里摘取一段网络代码：

```r
countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}

################################################################################
# An example
################################################################################
cnts <- c(4250, 3300, 200, 1750, 50, 0)
lens <- c(900, 1020, 2000, 770, 3000, 1777)
countDf <- data.frame(count = cnts, length = lens)

# assume a mean(FLD) = 203.7
countDf$effLength <- countDf$length - 203.7 + 1
countDf$tpm <- with(countDf, countToTpm(count, effLength))
countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
```

显而易见，要想通过 counts 计算 TPM，必须知道 count,efflength,但最难的也是这个有效长度

##### what is effective length of the feature?

I traced back to this paper by Lior pachter group [Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation](http://www.nature.com/nbt/journal/v28/n5/full/nbt.1621.html).

It is quite mathematical, but the [general idea is](https://groups.google.com/forum/#!searchin/kallisto-sleuth-users/Effective$20Length/kallisto-sleuth-users/SlJWXFMEEiM/ftkrtPZyAQAJ):

> If we take the fragment length to be fixed, then the effective length is how many fragments can occur in the transcript. This turns out to be length - frag_len +1. The number of fragments coming from a transcript will be proportional to this number, regardless of whether you sequenced one or both ends of the fragment. In turn, when it comes to probabilistically assigning reads to transcripts the effective length plays a similar role again. Thus for short transcripts, there can be quite a difference between two fragment lengths. To go back to your example if you have transcript of length 310, your effective length is 10 (if fragment length is 300) or 160 (if fragment length is 150) in either case, which explains the discrepancy you see.

From @Rob

> The effective length is computed by using the fragment length distribution to determine the effective number of positions that can be sampled on each transcript. You can think of this as convolving the fragment length distribution with the characteristic function (the function that simply takes a value of 1) over the transcript. For example if we observe fragments of length 50 --- 1000, a position more than 1000 bases from the end of the transcript will contribute a value of 1 to the effective length, while a position 150 bases will contribute a value of F(150), where F is the cumulative distribution function of the fragment length distribution. For single end data, where we can't learn an empirical FLD, we use a gaussian whose mean and standard deviation can be set with --fldMean and --fldSD respectively.

From Harold Pimentel's post above. He is in Lior Pachter's group.

> Effective length refers to the number of possible start sites a feature could have generated a fragment of that particular length. In practice, the effective length is usually computed as:

> where uFDL is the mean of the fragment length distribution which was learned from the aligned read. If the abundance estimation method you’re using incorporates sequence bias modeling (such as eXpress or Cufflinks), the bias is often incorporated into the effective length by making the feature shorter or longer depending on the effect of the bias.

所以说，自己计算 efflength 相当麻烦而且可能出错，有没有简单的方法呢？

有： 使用 salmon 软件会进行伪比对计算 TPM 和 count 同时会得到对应转录本的有效长度，将其计算结果中有效长度提取出来即可。

首先，提取对应有效长度给 fc 计算结果：

```python
import pandas as pd
import os,glob
fc_dir="lncRNA/normal/tiqu_trans"
sm_dir="LncRNA/salmon/normal/lncRNA"

fc=glob.glob(os.path.join(fc_dir,"*.lncRNA.trans.counts"))

for fs in fc:
    sample_name=fs.replace("lncRNA/normal/tiqu_trans/","").replace(".lncRNA.trans.counts","")
    df=pd.read_csv(fs,sep="\t",header=0,skiprows=1).iloc[:,[0,1]].set_index("Geneid")
    df.columns=['count']
    df2=pd.read_csv(os.path.join(sm_dir,"{}_quant".format(sample_name),"quant.sf"),sep="\t",header=0).set_index("Name")
    df1=df.join(df2[['EffectiveLength']])
    df1.to_csv(os.path.join(fc_dir,"{}.lncRNA.trans_eff.counts".format(sample_name)),sep="\t")

```

然后计算 TPM：

```r
options(scipen = 200)
options(digits = 3)
for (name in Sys.glob("*.lncRNA.trans_eff.counts")){
    countdata<-read.table(name,sep="\t",header = TRUE,row.names = "Geneid")
    countToTpm <- function(counts, effLen){
        rate <- log(counts) - log(effLen)
        denom <- log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }
    countdata<-na.omit(countdata)
    countdata$TPM<-with(countdata,countToTpm(countdata$count,countdata$EffectiveLength))
    id<-unlist(strsplit(name,"[.]"))[1]
    write.table(countdata,paste(id,".lncRNA.trans.TPM",sep = ""),quote = FALSE,sep="\t")
}

```

然后，将所有的 TPM 再次整合到表达矩阵中：

```python

import pandas as pd
import glob
df=pd.read_csv("sample.lncRNA.trans.TPM",sep="\t")
df=df[['TPM']]
for name in glob.glob("*.TPM"):
    sample=name.replace(".lncRNA.trans.TPM","")
    df2=pd.read_csv(name,sep="\t")
    df2=df2[['TPM']]
    df2.columns=[sample]
    df=df.join(df2)

df.drop(['TPM'],axis=1).to_csv("normal_lncRNA_trans_TPM.csv",sep="\t")
```

现在，又有一个问题了。……
基因的表达量怎么算呢，基因的有效长度没法算啊，因为 salmon 只能计算转录本的表达量啊。这里，我们可以通过转录本推算基因有效长度：

```python
import argparse
import pandas as pd
import numpy as np
import sys

def main(args):
    gtable = pd.read_csv(args.ginput,skiprows=1,sep="\t").set_index('Geneid')
    gtable = gtable.iloc[:,0:1]
    gtable.columns = ['Count']
    ttable = pd.read_table(args.tinput).set_index('Name')
    tgmap = pd.read_table(args.tgmap, names=['t', 'g']).set_index('t')
    gene_lengths = {}
    j = 0

    # Map over all gene groups (a gene and its associated transcripts)
    for g, txps in tgmap.groupby('g').groups.iteritems():
        if j % 500 == 1:
            print("Processed {} genes".format(j))
        j += 1
        # The set of transcripts present in our salmon index
        tset = []
        for t in txps:
            if t in ttable.index:
                tset.append(t)
        # If at least one of the transcripts was present
        if len(tset) > 0:
            # The denominator is the sum of all TPMs
            totlen = ttable.loc[tset,'TPM'].sum()
            # Turn the relative TPMs into a proper partition of unity
            if totlen > 0:
                tpm_fracs = ttable.loc[tset, 'TPM'].values / ttable.loc[tset,'TPM'].sum()
            else:
                tpm_fracs = np.ones(len(tset)) / float(len(tset))
            # Compute the gene's effective length as the abundance-weight
            # sum of the transcript lengths
            elen = 0.0
            for i,t in enumerate(tset):
                elen += tpm_fracs[i] * ttable.loc[t, 'EffectiveLength']
            gene_lengths[g] = elen

    # Give the table an effective length field
    gtable['EffectiveLength'] = gtable.apply(lambda r : gene_lengths[r.name] if r.name in gene_lengths else 1.0, axis=1)
    # Write it to the output file
    gtable.to_csv(args.output, sep='\t', index_label='Name')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute gene-level effective lengths")
    parser.add_argument('--ginput', type=str, help='gene level input table')
    parser.add_argument('--tinput', type=str, help='transcript level input table')
    parser.add_argument('--tgmap', type=str, help='transcript -> gene mapping')
    parser.add_argument('--output', type=str, help='output table with extra column')

    args = parser.parse_args()
    main(args)
```

```shell
python inferlength.py --ginput sample.lncRNA.counts --tinput sample.quant.sf --tgmap lnc_tx2gene.txt --output sample.eff.lnc.counts
```

然后，老办法计算 TPM,然后整合。

### 6.差异分析

差异分析有两种情况，一种是组间差异，一种是样本间差异，组间差异的意思就是，有两组样本做比较，每组样本默认其都是重复性实验，即两组进行比较；而样本间差异意思就是虽然也是两个组进行比较，但是每个组内的样本是独立的，之间也是有差异的。

对于两种情况，分别推荐组间差异分析使用**deseq2**;样本间差异分析使用**degseq**

通常情况下就是，使用 deseq2 捕获不到的差异，会被 degseq 得到，但是否符合心中预期，还需要进一步衡量，一般要做热图的话，degseq 得到的差异基因很难形成红白对照，聚类的话，样本和组很难聚在一起。

##### Deseq2:

```r
library(DESeq2)
path=getwd()
normal<-read.table(paste(path,"/","normal_lncRNA_smcounts.csv",sep = ""),header = TRUE,row.names = "Name",check.names=FALSE)
sick<-read.table(paste(path,"/","sick_lncRNA_smcounts.csv",sep = ""),header = TRUE,row.names = "Name",check.names=FALSE)
mydata<-cbind(normal,sick)
sample<-colnames(mydata)
group_lst=c(rep("normal",length(colnames(normal))),rep("sick",length(colnames(sick))))
colData <- data.frame(row.names=sample,group_list=group_lst)
dds <- DESeqDataSetFromMatrix(countData = round(mydata),colData = colData,design = ~ group_list)
register(MulticoreParam(4))
dds2<-DESeq(dds,parallel=TRUE)
saveRDS(dds2,paste(path,"/","lncRNA_normlized_by_deseq2_sm.rds",sep = ""))
#suppressMessages(dds2)
res <-  results(dds2, contrast=c("group_list","normal","sick"))
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
diff_gene_deseq <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.table(diff_gene_deseq,file = paste(path,"/","lncRNA_top_diff2_sm.txt",sep = ""),quote = FALSE,sep="\t")
plotdata<-counts(dds2, normalized = TRUE)
write.table(plotdata,file = paste(path,"/","lncRNA_normlized_by_deseq2_sm.txt",sep = ""),quote = FALSE,sep="\t")
```

### Degseq:

```r
library(qvalue)
library(samr)
library(impute)
library(matrixStats)
library("DEGseq")
path=getwd()
normal<-read.table(paste(path,"/","normal_lncRNA_smcounts.csv",sep = ""),header = TRUE,row.names = "Name",check.names=FALSE)
sick<-read.table(paste(path,"/","sick_lncRNA_smcounts.csv",sep = ""),header = TRUE,row.names = "Name",check.names=FALSE)
mydata<-cbind(normal,sick)
write.table(mydata,paste(path,"/","lncRNA_all_smcounts.txt",sep = ""),quote=FALSE,sep="\t")
matrix1 <- readGeneExp(paste(path,"/","lncRNA_all_smcounts.txt",sep = ""), geneCol=1, valCol=c(2:83))
matrix2 <- readGeneExp(paste(path,"/","lncRNA_all_smcounts.txt",sep = ""), geneCol=1, valCol=c(84:207))
DEGexp(geneExpMatrix1=matrix1, geneCol1=1, expCol1=c(2:83), groupLabel1="nomal",
               geneExpMatrix2=matrix2, geneCol2=1, expCol2=c(2:125), groupLabel2="sick",
                pValue=1e-3, zScore=4, qValue=1e-3, foldChange=2, thresholdKind=5,
               method="MARS", outputDir=path)
```

##### 合并 case，control 的 TPM 到一个矩阵：

```r
path=getwd()
normal<-read.table(paste(path,"/","normal_mRNA_smTPM.csv",sep = ""),row.names="Name",header = TRUE,check.names=FALSE)
sick<-read.table(paste(path,"/","sick_mRNA_smTPM.csv",sep = ""),row.names="Name",header = TRUE,check.names=FALSE)
mydata<-cbind(normal,sick)

write.table(mydata,"mRNA.salmon.TPM.txt",sep="\t",quote=FALSE)
```

##### 将表达量为 0 的过滤掉：

```r
path=getwd()
mydata<-read.table(paste(path,"mRNA.salmon.TPM.txt",sep="/"),header=TRUE)
#这里也可以在read.table里面添加check.names=FALSE
colnames(mydata)<-gsub("X","",colnames(mydata))
top_m<-read.table(paste(path,"top_mRNA_degseq_fliter.csv",sep="/"),header=TRUE,sep="\t")
top_m_TPM<-mydata[top_m$GeneNames,]
top_m_filter<-top_m_TPM[apply(top_m_TPM,1,mean)>0,]
write.table(top_m_filter,"top_m_filter.txt",sep="\t",quote=FALSE)
```

### 7.对获取到的差异转录本做 GO PATHWAY 富集分析

不推荐使用 Y 叔的`clusterprofile`,这个软件虽然易用，但是注释信息非常有限，功能也不完善，做富集分析和注释，还是需要使用在线的数据库注释，这里推荐几个非常好的网站：

- [DAVID](https://david.ncifcrf.gov/tools.jsp)
- [metascape](http://metascape.org/gp/index.html#/main/step1)
- [webgestalt](http://www.webgestalt.org/)
- [BioDbNet](https://biodbnet-abcc.ncifcrf.gov/)

**DAVID**这个网站被诟病许久，大家都在吐槽它更新缓慢的数据库，但是用起来还是真香，滑稽。

**metascape**是我用过的迄今为止最好用的网站，它可以接受几乎所有的 id，做的分析有一般分析和自定义分析，自定义分析里面可以帮助挑选出我们想要的通路，得到其对应的基因，这个真的非常好用，而且它有一个 excel 插件，可以直接在 Excel 里面运行，

**webgestalt**这个可能我还没发觉它真正的功能，我用着不好用，富集不全，但是图做的好，这个没得黑。

**BioDbNet**这个网站简直全能，可以做 id 转换，还可以做全注释，强无敌！！！

对 LncRNA 的注释，可以使用 R 包**LncPath**，这个包用法简单，但最好需要配合**WGCNA**使用

### 8.对得到的差异基因做可视化

目前我所做的可视化无非两种，一种是聚类热图，另一种是箱线图或者小提琴图，一种应用于大量基因，另一种是少量基因的情况。

热图很简单，直接用表达矩阵和分组信息做就可以了，默认参数，数据最好 log2(n+1)

箱线图同理，ggplot2 就可以做了

另外一种图是网络交互图，使用 igraph

### 9.寻找 lncRNA 的靶基因

针对 lncRNA 研究，其中比较关键的点是确定 lncRNA 的靶基因，如果能够有方法预测到 lncRNA 的靶基因，会减轻我们的工作量。而如何去测 lncRNA 的靶基因，可以从 lncRNA 不同的作用模式入手。

lncRNAs 作用范围广泛，机制非常复杂。根据 lncRNA 不同的作用模式比如顺式（cis）和反式（trans）之分，比如 Signal，Decoy，Duide，Scaffold，也可以根据 lncRNA 与不同的分子分为 DNA、RNA 和蛋白，总体上包括了转录和转录后水平。

##### 根据顺式（cis）和反式（trans）的作用模式来预测

RNA 聚合酶转录得到一个与调控蛋白相关的 LncRNA 从而影响附近区域的编码基因。

cis 作用靶基因预测，认为 LncRNA 的功能与其坐标临近的蛋白编码基因相关，位于编码蛋白上下游的 LncRNA 可能与启动子或者共表达基因的其他顺式作用元件有交集，从而在转录或者转录后水平对基因的表达进行调控。判断一个 LncRNA 具有 Cis 调控作用通常要同时满足以下几个条件：

- 附近的基因表达情况与其保持一致；
- 该基因失活后会影响周围基因的表达；
- 会影响附近同一位点的基因表达。

对于满足以上条件的 LncRNA，首先找出位于其上游或者下游附近（10K）的编码蛋白基因，通过对编码蛋白的功能富集分析，从而预测 LncRNA 的主要功能，为后续 Cis 作用分析打下基础。

LncRNA 与 DNA 结合蛋白相关联，并调控相关靶基因的表达。

Trans 作用靶基因预测基本原理认为 LncRNA 的功能跟编码基因的位置关系没有关系，而与其共表达的蛋白编码基因相关。也就是说，当 LncRNA 与一些距离较远的基因在表达量上存在正相关或负相关的情况时，就可以通过样本间 lncRNA 与蛋白编码基因的表达量相关性分析或 WGCNA 共表达分析来预测其靶基因。

需要注意的是，lncRNA 的 Trans 作用靶基因预测只适合于样本量大的情况，如果样本量太少（如 6 个以下），分析将不可靠。

所以，我的分析方法是，首先分析相关性，使用 p 值（FDR）和相关性参数确定相关性较强的基因作为备选

然后，分析 mRNA 上游 10K 和下游 20K 的 lncRNA 作为 cis；不在范围内的使用 RNAplex 计算结合能，结合能小于-30 作为 trans

**计算相关性：**

```r
library("fdrtool")
base<-getwd()
mydata<-read.table(paste(base,"top_m_filter.txt",sep="/"),header=TRUE,check.names=FALSE)
mydata<-as.data.frame(t(mydata))
chayi<-read.table(paste(base,"top_lnc_filter.txt",sep="/"),header=TRUE,check.names=FALSE)
chayi<-as.data.frame(t(chayi))
core<-function(x,y){
  cor.test(x,y)$estimate
}
pval<-function(x,y){
  cor.test(x,y)$p.value
}
fdr<-function(x){
  fdr<-fdrtool(x,statistic="pvalue")
  fdr$qval
}
mycor<-data.frame(test=rep(0,length(colnames(mydata))))
for (i in 1:length(colnames(chayi))){
  haha<-0
  for (j in 1:length(colnames(mydata))){
    xiangguan<-core(chayi[,i],mydata[,j])
    haha<-c(haha,xiangguan)
  }
  mycor<-cbind(mycor,as.data.frame(haha)[-1,])
}
mycor<-mycor[,-1]
colnames(mycor)<-colnames(chayi)
rownames(mycor)<-colnames(mydata)
write.table(mycor,paste(base,"correlation2.txt",sep="/"),quote=FALSE)
myp<-data.frame(test=rep(0,length(colnames(mydata))))
for (i in 1:length(colnames(chayi))){
  haha<-0
  for (j in 1:length(colnames(mydata))){
    xiangguan<-pval(chayi[,i],mydata[,j])
    haha<-c(haha,xiangguan)
  }
  myp<-cbind(myp,as.data.frame(haha)[-1,])
}
myp<-myp[,-1]
colnames(myp)<-colnames(chayi)
rownames(myp)<-colnames(mydata)
mfdr<-apply(myp,2,fdr)
write.table(mfdr,paste(base,"p_values(FDR)2.txt",sep="/"),quote=FALSE)

##guolv
haha<-data.frame("correlation"=0,"mRNA"=0,"lncRNA"=0)
xiangguan<-read.table(paste(base,"correlation2.txt",sep="/"),header = TRUE)
for (i in colnames(xiangguan)){
  if (length(row.names(subset(xiangguan,abs(xiangguan[,i])>0.6,select = i)))==0){
    next
  }else{
  xiang<-subset(xiangguan,abs(xiangguan[,i])>0.6,select = i)
  colnames(xiang)<-"correlation"
  xiang$mRNA<-rownames(xiang)
  xiang$lncRNA<-i
  rownames(xiang)<-NULL
  haha<-rbind(haha,xiang)
  }
}

heihei<-data.frame('p_values(FDR)'=0,"mRNA"=0,"lncRNA"=0)
colnames(heihei)[1]<-"p_values(FDR)"
xiangguan<-read.table(paste(base,"p_values(FDR)2.txt",sep="/"),header = TRUE)
for (i in colnames(xiangguan)){
  if (length(row.names(subset(xiangguan,xiangguan[,i]<0.05,select = i)))==0){
    next
  }else{
    xiang<-subset(xiangguan,xiangguan[,i]<0.05,select = i)
    colnames(xiang)<-"p_values(FDR)"
    xiang$mRNA<-rownames(xiang)
    xiang$lncRNA<-i
    rownames(xiang)<-NULL
    heihei<-rbind(heihei,xiang)
  }
}
haha<-haha[-1,]
heihei<-heihei[-1,]
fin<-merge(haha,heihei,by.x = c("lncRNA","mRNA"),by.y = c("lncRNA","mRNA"))
write.table(fin,paste(base,"xiangguanxing2.txt",sep="/"),quote = FALSE,row.names = FALSE)
```

> note: 这里计算的方法是使用 for 循环迭代结算，效率较低，不适合与大量数据，而且这里两个数据表一定要有相同的样本顺序，不然计算的相关性和 P 值会有显著的不同。最后，对于 FDR 矫正，这里只是针对一个 lncRNA 和多个 mRNA 的 P 值做了矫正。

```r
library(Hmisc)
library(reshape2)
library(fdrtool)
lncRNA<-read.table("lncRNA_difference.txt",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE,sep="\t")
mRNA<-read.table("mRNA_difference.txt",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE,sep="\t")
n_lnc<-nrow(lncRNA)
n_mRNA<-nrow(mRNA)
alldata<-rbind(lncRNA,mRNA)
mydata<-as.matrix(t(alldata))
p<-rcorr(mydata)
cor<-p[["r"]][1:n_lnc,(n_lnc+1):nrow(alldata)]
pvalue<-p[["P"]][1:n_lnc,(n_lnc+1):nrow(alldata)]
cor<-cbind(rownames(cor),cor)
colnames(cor)[1]<-"LncRNA"
pvalue<-cbind(rownames(pvalue),pvalue)
colnames(pvalue)[1]<-"LncRNA"
new_cor<-melt(as.data.frame(cor),id.vars="LncRNA",variable.name = "mRNA",value.name = "correlation")
new_p<-melt(as.data.frame(pvalue),id.vars="LncRNA",variable.name = "mRNA",value.name = "p.value")
jiaozheng<-function(x){fdrtool(x,statistic = "pvalue",plot = FALSE,verbose = FALSE)$qval}
all_cor<-merge(new_cor, new_p, by.x=c("LncRNA","mRNA"), by.y=c("LncRNA","mRNA"))
all_cor<-na.omit(all_cor)
all_cor$p.value<-as.numeric(all_cor$p.value)
all_cor$correlation<-as.numeric(all_cor$correlation)
all_cor$p.asjust<-jiaozheng(all_cor$p.value)
all_fli<-subset(all_cor,abs(all_cor$correlation)>0.6&all_cor$p.asjust<0.05)
write.table(all_fli,"xiangguanxing.txt",quote=FALSE,row.names=FALSE)
```

> 这个脚本可以解决迭代问题，而且避免了样本不一致的情况。
> **计算 cis**

首先，获取 bed 文件

```r
anno_m<-read.table("LncRNA/readscounts/gene_hit/hg19_all_transcript.bed",sep = "\t",header = FALSE)
anno_lnc<-read.table("LncRNA/readscounts/gene_hit/NONCODE.bed",sep = "\t",header = FALSE)
chayi<-read.table("xiangguanxing2.txt",header = TRUE)
lncRNA<-subset(anno_lnc,V4%in%chayi$lncRNA)
mRNA<-subset(anno_m,V4%in%chayi$mRNA)
write.table(lncRNA,"lncRNA.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(mRNA,"mRNA.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
```

所有基因的 bed 文件 noncode.bed 和 hg19.bed 可以在网上下载，也可以自己做，自己做的话使用 python 的 gtf_parse 包对 gtf 文件进行解析，然后做一个就行了

然后，使用 bedtools 计算位置信息：

```bash
windowBed -a mRNA.bed -b LncRNA.bed -l 10240 -r 20480 -sm > cis.txt
```

**计算 trans 结合能，需要首先获取转录组序列，然后运行 RNAplex**

```python
import pandas as pd
import subprocess,os,re
from Bio import SeqIO
hg19="LncRNA/testdir/all_mRNA_know.fa"
noncode="LncRNA/testdir/all_lncRNA_know.fa"
base="LncRNA/readscounts/salmon/gene_hit"
hg19_fasta={seq_record.id:str(seq_record.seq) for seq_record in SeqIO.parse(hg19,'fasta')}
lnc_fasta={seq_record.id:str(seq_record.seq) for seq_record in SeqIO.parse(noncode,'fasta')}



def get_energy(lncRNA,mRNA):


    fa="lnc_tmp.fasta"
    fa2="mRNA_tmp.fasta"
    with open(fa2,'w') as f:
        f.write(">{}\n{}".format(mRNA,hg19_fasta[mRNA]))
    with open(fa,'w') as l:
        l.write(">{}\n{}".format(lncRNA,lnc_fasta[lncRNA]))


    get_energy=['app/RNAplex-0.2/RNAplex/bin/RNAplex','-t',fa2,'-q',fa]

    cm3=subprocess.Popen(get_energy,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    out,err=cm3.communicate()
    ans=0
    for line in out.decode("utf-8").splitlines():
        if not line.startswith(">"):
            ans=re.findall("\(-?\w+.?\w+\)",line)[0].replace("(","").replace(")","")
    #print(ans)

    os.remove(fa2)
    os.remove(fa)
    return ans

with open("xiangguanxing2.txt") as f:
    for line2 in f:
        if line2.startswith("NONH"):
            lncRNA,mRNA,*ot=line2.split(" ")
            print("{}--{}:".format(lncRNA,mRNA),get_energy(lncRNA,mRNA))
```

> 注意：这里想要获取转录本序列，千万不要使用 bedtools getfast，因为这样获取的是整个转录本的基因组序列，包含内含子，所以，我们使用 cufflinks 软件中 gffread 命令来获取转录本序列

```bash
gffread  -w out.fa -g hg19.fa in.gtf
```

##### 根据 lncRNA 结合的核酸序列来预测

实际上这种方式是最容易想到的，因为 lncRNA 是核酸，而核酸与核酸之间是有碱基互补配对的，不管是 lncRNA 与 RNA（比如 microRNA，mRNA）之间还是 lncRNA 与 DNA 之间，我们可以根据互补配对这一特性进行预测。这里我们展开举三个例子：

**A、lncRNA-microRNA-mRNA**

microRNA 对 mRNA 靶基因预测的方式可以反过来用于 lncRNA 靶基因预测，这里适用的模式就是 lncRNA 作为 microRNA sponge 吸附 microRNA，或者进一步通过 ceRNA 的作用机制调控靶基因，原理都是 microRNA 与 lncRNA 和 mRNA 的序列结合，这种方式也是现在预测靶基因用的最多的，比如网站：starbase（http://starbase.sysu.edu.cn/）。

**B、lncRNA-mRNA**

在上面这种模式里，lncRNA 对靶基因 mRNA 的调控是通过 microRNA 来介导的，当然 lncRNA 与 mRNA 之间的直接互补也能帮我们预测靶基因，这种模式我们以前说过，比较常见的是反义 lncRNA（Antisense lncRNA）通过与 mRNA 结合形成 RNA 二聚体https://www.jianshu.com/p/92451bd7c030，保护 mRNA 使之更难被 RNA 酶降解。因此这种模式来预测靶基因的原理是预测 lncRNA 序列与 mRNA 序列结合的自由能，自由能越小，越容易结合。这种模式的在线工具有 lncRNATargets（http://www.herbbol.org:8001/lrt/）

**C、lncRNA-DNA**

在以前文章中我们也介绍过 lncRNA 通过结合单链 DNA 发挥调控作用的：

因此，基于同样的原理，理论上我们也可以进行预测，这种预测工具 lncRNATargets 同样可实现。

##### 根据蛋白的结合特性来预测:

lncRNA 核酸与核酸结合依据的是碱基配对原理，同样某些蛋白与核酸的结合也有一定规律，比如 RNA 结合蛋白，因此，我们可以依据这些规律预测 lncRNA 的结合蛋白，比如工具包括 RBPDB（http://rbpdb.ccbr.utoronto.ca/index.php）、RNAcommender（http://rnacommender.disi.unitn.it/）等。