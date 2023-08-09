---
title: 用VueJS和Fastapi通过websocket实现进度条追踪
date: 9th Aug 2023
description: 用VueJS和Fastapi通过websocket实现进度条追踪
image: https://cdn.sanity.io/images/8edntncj/production/7cb22a8a9885468d96662f6ece61ceb2bc95666d-500x300.png
alt: 用VueJS和Fastapi通过websocket实现进度条追踪
ogImage: https://cdn.sanity.io/images/8edntncj/production/7cb22a8a9885468d96662f6ece61ceb2bc95666d-500x300.png
tags: ['Python','前端','后端']
published: true
---

### 简介

在我们的应用小程序中，我们是前后端分离的。前端页面只负责渲染，而后端需要处理数据。但是如果遇到数据量很大的情况下，我们处理起来就很缓慢，如果我们想通过AJAX的方法追踪后台数据变化的进度，需要用到轮询的方案，这个是非常消耗资源的。这里我们用VueJS和Fastapi的小例子演示前端传递数据，后台用10秒处理数据并实时反应进度给前台的实现。

### 后端Fastapi

因为后端写起来相对简单，我们先写个后端

新建`my-fastapi`目录，然后新建`main.py`文件，代码如下所示：

```python
from time import sleep
from fastapi import FastAPI, WebSocket ##导入websocket
##创建app
app = FastAPI()
##路由端点需要websocket
@app.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()  ## 等待链接
    while True:
        data = await websocket.receive_json() ###等待接收Json数据
        ### 这是一个10S的任务
        for i in range(10): 
            sleep(1)
            ### 将此刻的循环次数传递给前端
            await websocket.send_json({"data": "", "time": (i+1)*10}) 
        jieguo = data['input'] + " finished!"
        ### 等待返回数据
        await websocket.send_json({"data": jieguo, "time": 100})
```

### 前端

前端用Vuejs+NaiveUI实现：

创建一个Vite项目：

```bash
npm create vite@latest my-vue-app -- --template vue
```

安装 Naive-UI

```bash
npm i -D naive-ui vfonts
```

将`src/app.vue` 中的内容删除，替换成下面的：

```vue
<script setup>
///UI
import { NButton, NInput, NSpace, NH1, NProgress, useThemeVars } from 'naive-ui'
import { onMounted, ref } from 'vue'
import { changeColor } from "seemly";
/// 绑定数据
const data = ref('')
const tm = ref('')
const inputData = ref('')
/// 建立持久通信
const connection = new WebSocket("ws://localhost:8000/ws")
const themeVars = useThemeVars()

///这里用来向后端发送Json数据，注意要转换为字符串格式 JSON.stringify()
const submit = () => {
  connection.send(JSON.stringify({ input: inputData.value }))
}
/// 这里接收数据，注意要把字符串格式的JSON解析为json对象 JSON.parse()
onMounted(() => {
  connection.onmessage = function (e) {
    const backData = JSON.parse(e.data)
    data.value = backData.data
    tm.value = backData.time
  }

})
</script>

<template>
  <n-space justify="center" vertical>
    <n-h1 style="text-align: center;">Hello {{ data }}</n-h1>
    <n-space  justify="center">
      <n-progress style="width: 300px;"  type="line" indicator-placement="inside" :color="themeVars.errorColor"
        :rail-color="changeColor(themeVars.errorColor, { alpha: 0.2 })" :percentage="tm" 
        />
    </n-space>
    <n-space style="display: flex; margin: 10px auto; justify-content: center;">
      <n-input v-model:value="inputData" type="text" placeholder="测试输入" />
      <n-button type="primary" @click="submit">Submit</n-button>
    </n-space>
  </n-space>
</template>
```

最终的效果如图所示：

*这个动图没有录制完整，最终会到100%。*

![showing](https://picshack.net/ib/Psoo0yimK7.gif)