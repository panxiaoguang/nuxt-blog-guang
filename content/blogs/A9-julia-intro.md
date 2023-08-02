---
title: Julia 入门指难
date: 8th May 2021
description: Julia 入门指难.
image: /blogs-img/julia.webp
alt: Julia 入门指难
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

### 如何安装Julia

有很多方法，其中最简单的就是去各大景象站点下载编译好的二进制包，例如
- [清华大学开源软件镜像站](https://mirrors.tuna.tsinghua.edu.cn/help/julia-releases/)
- [北京外国语大学开源软件镜像站](https://mirrors.bfsu.edu.cn/help/julia-releases/)
- [上海交通大学软件源镜像服务](https://mirrors.sjtug.sjtu.edu.cn/julia-releases/)

另外，可以使用包管理工具[jill](https://github.com/johnnychen94/jill.py)下载安装，

- 安装/更新 jill： `pip install jill --user -U` (需要 Python 3.6 或更新的版本)

- 安装 Julia：`jill install [VERSION] [--upstream UPSTREAM] [--confirm]`

- 查询现存的上游镜像：`jill upstream`

- 帮助文档：`jill [COMMAND] --help`

利用 jill 安装完成后即可通过在命令行执行 `julia/julia-1/julia-1.4` 来启动不同版本的 Julia

当然，你也可以自己从源码构建安装，但是会很慢，并不推荐。

### 如何安装包

julia 包管理器类似于R，可以进入REEL中，然后输入`]`，此时进入了安装包模式，输入`add {packages}`即可完成安装

### 安装包的速度很慢？

因为Julia的包托管在github上，慢是很正常的，此时你需要一份镜像代理。

具体使用方法可以参考[清华镜像julia使用指南](https://mirrors.tuna.tsinghua.edu.cn/help/julia/)

## 如何查找有哪些包？

我列举几个网站，可以从这些站点寻找已经注册的包，小白用户就不要考虑未注册的包了。

- [julia-hub](https://juliahub.com/)
- [julia-observer](https://juliaobserver.com/)
- [julia-packages](https://juliapackages.com/)

### 有没有比较系统的教程？

有，[网易云课堂](https://study.163.com/course/introduction.htm?courseId=1208959805&_trace_c_p_k2_=d50856be95c5489c94f5fa8f276bc000)
有免费的教学视频，相对全面且新手友好！版本也较新

### julia可以使用jupyter notebook 吗？

可以，在julia中安装`IJulia`包即可

### 可以买书学习吗？

不推荐，因为julia在1.0版本前和版本后就像两个语言，有很大的变化，买的书很可能不起作用反而给你带来困扰。即使要买，也要看好版本，尽量选择1.0版本后的。

### julia可以用来做生信分析吗？

完全可以，`BioJulia`了解一下! 但是目前转录组的相关包还没有出现。
转录组还是用R包吧！

### 有问题去哪里问？

julia有官方中文论坛，可以去那里提问问题，地址是
[julia中文社区](https://discourse.juliacn.com/)
你也可以加入官方QQ群：316628299

### 什么情况下我要学习Julia

当你再也无法忍受python的蜗牛速度，而且无法通过算法优化来提升性能，你可以果断尝试Julia，或许会得到意想不到的速度提升。