---
title: python, perl 和julia的性能对比
date: 7th May 2021
description: python, perl 和julia的性能对比.
image: /blogs-img/language.webp
alt: python, perl 和julia的性能对比
ogImage: /blog-img/language.webp
tags: ['r语言', 'python', 'julia']
published: true
---

*2023/3/20更新*：
>Codon是一个高性能的Python编译器，它将Python代码编译为本地机器代码，而不需要任何运行时开销。Python上的典型加速在单个线程上大约为10-100x或更多。Codon的性能通常与C/C++不相上下。与Python不同，Codon支持本机多线程，这会导致速度提高很多倍。Codon可通过插件基础设施进行扩展，使您能够合并新的库、编译器优化甚至关键字。

现在，让我们测试codon是否能给python提速，在此之前，我们需要修改以下python的代码

```python
import sys
def calculateGC(sequence:str)->Tuple[int,int]:
    """Calculate the GC content of a DNA sequence"""
    gc = 0
    allnumber = 0
    for i in sequence:
        if i != 'N' and i != 'n':
            allnumber += 1
            if i == 'G' or i == 'C' or i == 'g' or i == 'c':
                gc += 1
    return gc, allnumber


def main(file:str):
    """Main function"""
    gcNum = 0
    allNum = 0
    with open(file,'r') as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                continue
            else:
                gc, allnumber = calculateGC(line)
                gcNum += gc
                allNum += allnumber
    print(f'GC content is: {gcNum/allNum:.3f}')

main(sys.argv[1])
```
然后运行
```
codon build --release -o calGC calGC.py
```
最后速度为：

```bash
time ./calGC ../../Project/DataBase/hg38.fa     
GC content is: 0.410
./calGC ../../Project/DataBase/hg38.fa  22.30s user 2.75s system 111% cpu 22.421 total
```
速度还不错，已经可以超越不用Biojulia包的julia函数了。

*2022/11/14更新*：

最近python3.11出来了，据说性能有很大的提升, 一位国外的小哥(Dennis Bakhuis)采用简单的蒙特卡洛预测圆周率的方式测试循环的性能，发现Python3.11性能确实是突飞猛进，同一时间，Julia也已经更新到V1.8了，于是我在他的github下贡献了julia版本的代码，希望继续比较多个语言的计算特性。

总而言之，虽然python性能进一步优化，但和julia相比，速度依旧不够打，1000000个循环，python3.11用时6秒，而julia仅需要0.033秒。具体可以看github
[https://github.com/dennisbakhuis/python3.11_speedtest](https://github.com/dennisbakhuis/python3.11_speedtest)

> 在生物信息学中经常用到的脚本语言主要是`python`和`perl`，他们被用来**处理文本**，**大量统计**，**流程控制**等等，其自身也是各有优势。比如说`perl`天生就为了处理文本而生，但是`python`确是有名的胶水语言，特别在整合`C`代码时显示出巨大的优势，其语法简洁易懂，易于维护更让其成为仅次于`C`和`JAVA`的第三大语言，但其糟糕的性能在处理大量循环时会让人忍不住抓狂。因此，`Julia`语言应运而生，其控制了`python`中没必要的动态性，加之使用JIT技术让其能够保有高性能的同时具备简洁的语法。

>说了那么多，在生物信息上我们经常需要处理大量的文本文件，例如`Fasta`格式的序列文件，那么三者又是谁快呢？

### 版本控制
- python3 = 3.8.3
- perl = 5.26.2
- julia = 1.5.0-beta
- system = centos 8

### 计算内容
从UCSC上下载人类参考基因组 [hg38.fa.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) 并解压，计算基因组GC含量，**N碱基不算在总长中**。

### 代码
**perl**
```perl
#!/usr/bin/perl -w

use strict;

if(@ARGV < 1){
    die "Usage : perl $0 <genome.fa>\n";
}

my $input = shift @ARGV;

my ($sum,$G_num,$C_num,$N_num)=(0,0,0,0);
my $id;
open IN, "< $input" or die $!;
while(my $line = <IN>){
    chomp $line;
    if($line =~ />([^\s]+)/){
        $id = $1;
    }else{
        $sum += length($line);
        $G_num += ($line =~ tr/Gg/Gg/);
        $C_num += ($line =~ tr/Cc/Cc/);
        $N_num += ($line =~ tr/Nn/Nn/);
    }
}
close IN;

my $GC_rate = ($G_num+$C_num)/($sum-$N_num);
printf "GC content: %.3f \n",$GC_rate;
```

**julia**

```julia
function lineGC(seq::String)
    GCnumber=count(x->(x=='G'||x=='C'||x=='g'||x=='c'),seq)
    lineNum=count(x->(x!='N' && x!='n'),seq)
    (GCnumber,lineNum)
end

function calGC(fs)
    GCnumber=zero(Int)
    lineNum=zero(Int)
    open(fs,"r") do IOstream
        for line in eachline(IOstream)
            if startswith(line,">")
                continue
            else
                GC,all=lineGC(line)
                GCnumber+=GC
                lineNum+=all
            end
        end
    end
    round(GCnumber/lineNum;digits=3)
end

println("GC content: ",calGC(ARGS[1]))
```

**python**

```python
import sys

def lineGC(seq):
   tmp=[base for base in seq if base =="G" or base =="g" or base == "C" or base == "c"]
   gcNumber=len(tmp)
   tmp2=[base for base in seq if base !="N" and base !="n"]
   allNumber=len(tmp2)
   return (gcNumber,allNumber)


with open(sys.argv[1],'r') as f:
    gcNum=0
    allNum=0
    for line in f:
       if line.startswith(">"):
           continue
       else:
           gc,alln=lineGC(line.strip("\n"))
           gcNum=gcNum+gc
           allNum=allNum+alln

print("GC content: {:.3f}".format(gcNum/allNum))
```

### 运行时间测试

python

![python](https://s1.ax1x.com/2020/09/05/wEU7ut.md.png)

julia

![julia](https://s1.ax1x.com/2020/09/05/wEae29.md.png)

perl

![perl](https://s1.ax1x.com/2020/09/05/wEatxA.md.png)

### 总结

结果令人咋舌，可以从sys时间看出来python和perl都是立马启动，而julia在函数的即时编译上花了一点时间（一半时间）。 总体用时上，julia仅比perl快了1秒，而python却用了惊人的9分钟，😭
### 后记

python 也不是这么不堪，想要提速还是可以有很多办法的，比如切换pypy, 或者也用正则表达式，例如：

```python
import sys
import re


def lineGC(seq):
    pattern_1 = re.compile(r"G|C",re.I)
    pattern_2 = re.compile(r"N",re.I)
    gcNumber=len(pattern_1.findall(seq))
    allNumber=len(seq)-len(pattern_2.findall(seq))
    return (gcNumber,allNumber)


with open(sys.argv[1],'r') as f:
    gcNum=0
    allNum=0
    for line in f:
       if line.startswith(">"):
           continue
       else:
           gc,alln=lineGC(line.strip("\n"))
           gcNum=gcNum+gc
           allNum=allNum+alln

print("GC content: {:.3f}".format(gcNum/allNum))
```

> 这样计算下来，大概需要6分20秒，提速了一半

另外，我们也可以使用NumPy的向量化运算来提速

```python
import sys
from pyfaidx import Fasta
import numpy as np

def lineGC(seq):
    gc_number = np.where((seq==b'G')|(seq==b'C')|(seq==b'g')|(seq==b'c'))[0].shape[0]
    n_number = np.where((seq==b'N')|(seq==b'n'))[0].shape[0]
    allnumber = seq.shape[0] - n_number
    return (gc_number,allnumber)

def calGC(fs):
    GC = 0
    all = 0
    hg38 = Fasta(fs)
    for record in hg38:
        seq = np.asarray(record)
        gc_number,all_number=lineGC(seq)
        GC = GC + gc_number
        all = all + all_number
    return (GC, all)

if __name__ == "__main__":
    gcNum, allNum = calGC(sys.argv[1])
    print("GC content: {:.3f}".format(gcNum/allNum))
```
> 这样的话，就只需要2分22秒了，已经是非常快的了，但是和perl还是有差距的。

**最后，难道julia真的速度和perl就相差无几吗?**

答案是否定的，因为julia设计是为了科学计算的，但是其字符串的性能并算不上优秀，我们可以调用`BioSequence`来处理生物序列

```julia
using BioSequences
using FASTX

function lineGC(seq)
    GCnumber=count(x->(x==DNA_G||x==DNA_C),seq)
    lineNum=length(seq)-count(isambiguous,seq)
    GCnumber,lineNum
end

function calGC(fs)
    GCnumber=zero(Int)
    lineNum=zero(Int)
    reader=open(FASTA.Reader,fs)
    for record in reader
        GC,all=lineGC(FASTA.sequence(record))
        GCnumber+=GC
        lineNum+=all
    end
    close(reader)
    round(GCnumber/lineNum;digits=3)
end

println("GC content: ",calGC(ARGS[1]))
```

> 这样就只需要**11**秒就可以计算出答案了