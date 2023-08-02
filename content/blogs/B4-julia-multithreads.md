---
title: julia 多线程
date: 8th May 2021
description: julia 多线程.
image: /blogs-img/julia.webp
alt: julia 多线程
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

`julia`本身是一门很快速的语言，但是现代计算机往往具有多核心多线程设计，因此，充分发挥硬件，能进一步提高效率

### 多线程的启动

`julia`从1.5开始新添加了命令行参数`-t num_procs`，例如，你想启动一个10个线程的julia，那么就可以执行：

```bash
julia -t 10
```

进入`REPL`后，可以查看当前线程数：
```julia
Threads.nthreads()
10
```
### 如何多线程？

`julia`提供了一个简单的宏函数，来快速的包装`for`循环，使其可以将任务拆分为`n`份，每分独立执行


```julia
Threads.@threads
```

例如：


```julia
a = zeros(10)

Threads.@threads for i = 1:10
    a[i] = Threads.threadid()
end
```

非常简单的，可以快速将你的单线程函数转化为多线程。

### 线程安全

举个例子：

```julia
acc = Ref(0)

@threads for i in 1:1000
    acc[] += 1
end
```
正常情况下，该流程计算的是从累加1000个1，结果应该是1000才对，但是因为多线程调度，导致acc最终的结果并非为1000，即得到了错误的计算结果。

第一个解决办法是使用**原子操作**

```julia
acc = Atomic{Int64}(0)

@threads for i in 1:1000
    atomic_add!(acc, 1)
end
```

另一种方法是在变量上添加线程锁，使其同一时间有且只有一个线程可以操作该变量

```julia
acc = Ref(0)
l = ReentrantLock()

@threads for i=1:1000
    lock(l) do 
        acc[]+=1
    end
end
```


当然，举例只是为了说明多线程的安全问题，该变量的累积添加了线程锁后其实和单线程操作没有任何区别，反而因为线程调度影响了性能。因此，需要综合考虑你的项目是否需要多线程。