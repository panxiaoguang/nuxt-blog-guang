---
title: Julia 计算相关性检验
date: 7th Oct 2021
description: Julia 计算相关性检验.
image: /blogs-img/julia.webp
alt: Julia 计算相关性检验
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

众所周知，计算相关性非常的简单，因为`R` 语言中有函数`cor.test()`,该函数可以计算多种方法的相关性检验，返回相关性，Pvalue等检验值，但是这个函数在`Julia`中并不存在，让Julia作为一门科学计算语言显得并不完美。

首先，我们来看一下R语言计算线性相关的结果是什么样子的？

```R
a<-c(1,3,4,5,7)
b<-c(2,4,6,8,9)
c<-cor.test(a,b)
c
```

```R
	Pearson's product-moment correlation

data:  a and b
t = 7.7771, df = 3, p-value = 0.004423
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.6757774 0.9984873
sample estimates:
     cor 
0.976086
```

首先，我们生成了两个示例的向量，这两个向量必须`等长`

然后，我们使用了`cor.test`，因为默认就是皮尔森，因此我们不再使用method参数，得到的结果也很丰富。


但是，我们比较Julia呢？

首先，julia标准库并没有计算相关性的函数，这里我们可以调用Statistics包来计算相关性

```Julia
a = [1,3,4,5,7]
b = [2,4,6,8,9]

using Statistics

c = cor(a,b)

0.9760860118037878
```
然后，我们参考R 语言的定义计算p-value和95%的置信区间。

```Julia
using Distributions ##载入包，用来计算T统计和t分布的cdf值
using Statistics ##计算cor

## 根据r和样本量n计算置信区间
function interval(r::Float64, n::Int64)
    z = 0.5 * log(exp(1), (1 + r) / (1 - r))
    s = 1 / (sqrt(n - 3))
    zl = z - 1.96 * s
    zh = z + 1.96 * s
    low = (exp(2 * zl) - 1) / (exp(2 * zl) + 1)
    high = (exp(2 * zh) - 1) / (exp(2 * zh) + 1)
    (low, high)
end

## 输入两个向量，返回各个参数
function cor_test(x::AbstractVector, y::AbstractVector)
    if length(x) == length(y)
        r = cor(x, y) #相关性值
        n = length(x) 
        df = n - 2 ## 自由度
        t = r * sqrt((df) / (1 - r^2)) ##T统计量
        tdist = cdf(TDist(df), t) ##t分布的累计
        if r > 0
            pvalue = (1 - tdist) * 2
        else
            pvalue = (tdist) * 2  
        end
    else
        println("请输入等长的两个向量")
    end
    fin = (r = r, pvalue = pvalue, df = df, t = t, cfi = interval(r, n))
end
```
测试以下函数的返回是否和R语言一致

```Julia
cor_test(a,b)

(r = 0.9760860118037878, pvalue = 0.004423314933107214, df = 3, t = 7.777137710478182, cfi = (0.6757635765936858, 0.9984873283114327))

c = cor_test(a,b)

c.r
0.9760860118037878

c.pvalue
0.004423314933107214
```