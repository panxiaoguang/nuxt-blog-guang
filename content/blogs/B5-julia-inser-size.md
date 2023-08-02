---
title: 用julia语言计算测序数据的Insert Size?
date: 5th Feb 2023
description: 用julia语言计算测序数据的Insert Size.
image: /blogs-img/julia.webp
alt: 用julia语言计算测序数据的Insert Size?
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

### Julia读取BAM的库
想要计算`Insert size`，需要提供一个基因组比对后的文件，`sam`也好，`bam`也罢。那么，使用julia语言计算该值的第一步便是了解如何读取和解析`BAM`文件格式。


我们使用`BioJulia`提供的`XAM`包来读取BAM文件。所以我们需要首先安装该包。
打开julia，输入`]`进入Pkg模式
```r
add XAM
```

### 使用XAM读取和解析BAM文件的一般格式

```r
reader = open(BAM.Reader, "data.bam")
record = BAM.Record()
while !eof(reader)
    empty!(record)
    read!(reader, record)
    # 做一些事情，例如解析和计算
end
```

上面的流程总体上就是
- 使用BAM.Reader读取文件
- 定义一个空的Record对象
- 从头至尾循环整个BAMfile
- 将每行bam读入Record

这样的操作方式可以节省内存，避免循环很大的bam文件爆内存。

### 什么是Insert Size？

通俗的讲，Insert 长度就是指双端序列比对后，模版的长度。所以我们要计算需要保证如下条件

1. reads是成对，最好是Proper paired
2. 只需要计算read1即可，不然就算重了
3. 即使Proper paired也会有负的Insert，需要移除

代码如下

```r
using XAM
using Statistics

function insert_size_dist(Reader::XAM.BAM.Reader)
    insert_length = Int64[]
    record = BAM.Record()
    while !eof(Reader)
        empty!(record)
        read!(Reader, record)
        if BAM.flag(record) & 0x2 != 0 ## paired
            if BAM.flag(record) & 0x40 != 0  ## first in pair
                t_len = BAM.templength(record)
                if t_len > 0
                    push!(insert_length, t_len)
                end
            end
        end
    end
    (mean(insert_length), std(insert_length))
end
```

