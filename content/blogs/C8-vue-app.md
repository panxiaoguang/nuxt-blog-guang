---
title: 生信小白如何用vuejs简单创建一个交互式网页
date: 8th Aug 2023
description: 生信工作者常常需要创建交互式App，这通常是一个网站，但是对于非专业的前端来说，如何更简单完成这个任务呢？下面由我来用一个例子带你认识这个过程！
image: https://cdn.sanity.io/images/8edntncj/production/7cb22a8a9885468d96662f6ece61ceb2bc95666d-500x300.png
alt: 生信小白如何用vuejs简单创建一个交互式网页
ogImage: https://cdn.sanity.io/images/8edntncj/production/7cb22a8a9885468d96662f6ece61ceb2bc95666d-500x300.png
tags: ['前端']
published: true
---

### 什么是Vue.js

VueJS是一个渐进式的前端框架，所谓渐进式的意思就是你可以用它快速完成原型创作，然后在此基础上逐步完善。他可以足够简单，也可以足够完善，那么对于新手小白来说，这简直就是福利！:grin:

[Vuejs](https://cn.vuejs.org/) 相当于使用nodejs来渲染一个html模板，从而动态的生成我们想要的界面（DOM），类似的功能我们可能在Python中见过，例如Jinja2 或者 Django template. 就是给定一个基础的模板，然后我们只需要把给定的数据输入，就可以得到想要的页面了！

### 创建项目工程目录

这里我们使用[Vite](https://cn.vitejs.dev)来开发前端页面。创建一个Vue的模板

```bash
npm create vite@latest my-vue-app -- --template vue
```

安装[NaiveUI](https://www.naiveui.com/zh-CN/os-theme/docs/installation)

```bash
npm i -D naive-ui vfonts
```

配置NaiveUI

因为NaiveUI是可以任意引用的(Tree shaking),所以我们直接饮用想使用的code即可，无需配置。

### 创建单页面应用

在SRC目录中新建一个Components页面，或者Views 或者 Pages,什么名字都可以，我们把该页面的逻辑放进去。

例如，我们创建了`/Components/upload.vue`.

```vue
<script setup>
import { ref } from "vue";
import { lyla } from "lyla";
import { useMessage } from "naive-ui";
const message = useMessage();
const DataFromServer = ref([]);
const sampleName = ref(4);
const barcode = ref(23);
const chip = ref(27);
const lane = ref(28);
const datapath = ref(48);
const path = ref("10.2.100.1:/pakpox/pob8d1/");
const code1 = ref("");
const code2 = ref("");
const tableUpload = ref();
const params = ref({
  sampleName: sampleName,
  barcode: barcode,
  chip: chip,
  lane: lane,
  datapath: datapath,
  path: path,
});
const params2 = ref({
  file: "",
  sampleName: sampleName,
  barcode: barcode,
  chip: chip,
  lane: lane,
  datapath: datapath,
  path: path,
});
const customUpload = ({
  file,
  data,
  headers,
  withCredentials,
  action,
  onFinish,
  onError,
  onProgress,
}) => {
  const formData = new FormData();
  params2.value.file = file.file;
  if (data) {
    Object.keys(data).forEach((key) => {
      formData.append(key, data[key]);
    });
  }
  formData.append("file", file.file);
  lyla
    .post(action, {
      withCredentials,
      headers,
      body: formData,
      onUploadProgress: ({ percent }) => {
        onProgress({ percent: Math.ceil(percent) });
      },
    })
    .then(({ json }) => {
      if (json.hasBlank == 1) {
        message.warning("Upload success, but some blank in table");
      }
      if (json.hasDup == 1) {
        message.warning("Upload success, but some duplicate in table");
      }
      DataFromServer.value = json.rawData;
      code1.value = json.cmd1;
      code2.value = json.cmd2;
      onFinish();
    })
    .catch((error) => {
      message.success(error.message);
      onError();
    });
};
const customUpload2 = ({
  file,
  data,
  headers,
  withCredentials,
  action,
  onFinish,
  onError,
  onProgress,
}) => {
  const formData = new FormData();
  params2.value.file2 = file.file;
  if (data) {
    Object.keys(data).forEach((key) => {
      formData.append(key, data[key]);
    });
  }
  formData.append("file2", file.file);
  lyla
    .post(action, {
      withCredentials,
      headers,
      body: formData,
      onUploadProgress: ({ percent }) => {
        onProgress({ percent: Math.ceil(percent) });
      },
    })
    .then(({ json }) => {
      if (json.hasBlank == 1) {
        message.warning("Upload success, but some blank in table");
      }
      if (json.hasDup == 1) {
        message.warning("Upload success, but some duplicate in table");
      }
      DataFromServer.value = json.rawData;
      code1.value = json.cmd1;
      code2.value = json.cmd2;
      onFinish();
    })
    .catch((error) => {
      message.success(error.message);
      onError();
    });
};
const columns = [
  {
    title: "sampleName",
    key: "sampleName",
  },
  {
    title: "barcode",
    key: "barcode",
  },
  {
    title: "chip",
    key: "chip",
  },
  {
    title: "lane",
    key: "lane",
  },
  {
    title: "dataPath",
    key: "dataPath",
    width: 100,
    ellipsis: {
      tooltip: true,
    },
  },
];
const submitdata = () => {
  tableUpload.value?.submit();
};
</script>
<template>
  <n-layout>
    <n-layout-header bordered>
      <div style="margin: 10px 12px;">
      <span style="font-size: 1.3rem;font-weight: 500;">Glims to copy</span>
    </div>
    </n-layout-header>
    <n-layout-content>
      <div style="margin-top: 20px;">
      <n-grid x-gap="12" :cols="12">
        <n-gi :span="1"> </n-gi>
        <n-gi :span="3">
          <n-card hoverable>
            <n-space vertical>
              <n-p>Excel from Glims...</n-p>
              <n-upload
                ref="tableUpload"
                :default-upload="false"
                action="https://glims2excel-1-w4936186.deta.app/upload"
                response-type="json"
                :custom-request="customUpload"
                :data="params"
              >
                <n-button style="height: 35px">Upload..</n-button>
              </n-upload>

              <n-p>Index for sampleName</n-p>
              <n-input-number v-model:value="sampleName" clearable />
              <n-p>Index for barcode</n-p>
              <n-input-number v-model:value="barcode" clearable />
              <n-p>Index for chip</n-p>
              <n-input-number v-model:value="chip" clearable />
              <n-p>Index for lane</n-p>
              <n-input-number v-model:value="lane" clearable />
              <n-p>Index for dataPath</n-p>
              <n-input-number v-model:value="datapath" clearable />
              <n-p>Remoth path</n-p>
              <n-input v-model:value="path" type="text" />
              <n-button type="error" @click="submitdata"> Submit! </n-button>
              <n-divider dashed> Change Names </n-divider>
              <n-upload
                action="https://glims2excel-1-w4936186.deta.app/changenames"
                response-type="json"
                :custom-request="customUpload2"
                :data="params2"
              >
                <n-button round type="primary">Upload..</n-button>
              </n-upload>
            </n-space>
          </n-card>
        </n-gi>
        <n-gi :span="7">
          <n-space vertical style="margin-top: 50px">
            <n-tabs type="segment">
              <n-tab-pane name="table" tab="Table Preview">
                <n-data-table
                  :columns="columns"
                  :data="DataFromServer"
                  max-height="600px"
                />
              </n-tab-pane>
              <n-tab-pane name="cmd" tab="Shell cmd">
                <n-tabs type="line" animated>
                  <n-tab-pane name="shell1" tab="Shell1">
                    <n-card style="overflow-y: scroll; max-height: 700px">
                      <n-code :code="code1" language="bash" word-wrap> </n-code>
                    </n-card>
                  </n-tab-pane>
                  <n-tab-pane name="shell2" tab="Shell2">
                    <n-card style="overflow-y: scroll; max-height: 700px">
                      <n-code :code="code2" language="bash" word-wrap> </n-code>
                    </n-card>
                  </n-tab-pane>
                </n-tabs>
              </n-tab-pane>
            </n-tabs>
          </n-space>
        </n-gi>
        <n-gi :span="1"> </n-gi>
      </n-grid>
      </div>
    </n-layout-content>
    <n-layout-footer> 
      <div style="height: 90px;text-align: center;padding-top: 30px;">
        <p>built by XiaoguangPan</p>
        <p>Based on Vue.js & Naive UI</p>
      </div>
       </n-layout-footer>
  </n-layout>
</template>
```
然后我们把该组件引入App.vue,

```vue
<script setup>
import upload from "./components/upload.vue";
import hljs from 'highlight.js/lib/core'
import bash from 'highlight.js/lib/languages/bash'
hljs.registerLanguage('bash', bash)
</script>

<template>
  <n-message-provider>
    <n-config-provider :hljs="hljs">
      <upload />
    </n-config-provider>
    
  </n-message-provider>
</template>
```

就这么简单，短短几行代码，一个前端页面完成。完整的项目可以在[github](https://github.com/panxiaoguang/my_vue_learning)访问。

### 部署单页面应用到github page

因为我们没有后端，完全的就是一个客户端渲染页面，所以可以按照纯静态去部署。

在当前目录创建`.github/workflows/deploy.yml`文件

将以下流程代码写入：

```yaml
# 将静态内容部署到 GitHub Pages 的简易工作流程
name: Deploy static content to Pages

on:
  # 仅在推送到默认分支时运行。
  push:
    branches: ['main']

  # 这个选项可以使你手动在 Action tab 页面触发工作流
  workflow_dispatch:

# 设置 GITHUB_TOKEN 的权限，以允许部署到 GitHub Pages。
permissions:
  contents: read
  pages: write
  id-token: write

# 允许一个并发的部署
concurrency:
  group: 'pages'
  cancel-in-progress: true

jobs:
  # 单次部署的工作描述
  deploy:
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v3
      - name: Set up Node
        uses: actions/setup-node@v3
        with:
          node-version: 18
          cache: 'npm'
      - name: Install dependencies
        run: npm install
      - name: Build
        run: npm run build
      - name: Setup Pages
        uses: actions/configure-pages@v3
      - name: Upload artifact
        uses: actions/upload-pages-artifact@v1
        with:
          # Upload dist repository
          path: './dist'
      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v1
```

### 总结

我们借助Vuejs 和 NaiveUI, 完全不需要我们知道Html和CSS的知识，我们只需要把需要的组件调用即可。同时，我们的后端API已经在上一篇文章写好部署，我们直接调用即可。








