---
title: 把Fastapi部署到免费的Deta Space上
date: 6th Aug 2023
description: 我们去 “白嫖” Deta Space的功能来完成后端的部署，这样我们的前端在任何地方都可以调取该API，前端就是个客户端APP而已。
image: https://picshack.net/ib/v47sJyMyyy.png
alt: 把Fastapi部署到免费的Deta Space上
ogImage: https://picshack.net/ib/v47sJyMyyy.png
tags: ['Python','后端']
published: true
---

上一篇文章我们写了一个Streamlit的程序来全栈的执行我们的任务，但是我们也看到了它的一个缺点：**前端界面非异步，UI定制缺乏灵活性**。那么，我们接下来尝试采用**前后端分离**的方式来完成[上次](/blogs/C6-streamlit-data-app.md)的任务。

- 前端：Vue.js + Naive-UI + Vite + github page
- 后端：Fastapi

那么我们可以看到，我们是需要前端通过AJAX的方式从后端获取数据，后端异步的执行我们表格的处理任务，并生成脚本，推送给前端。
要想执行后端任务，必然需要服务器。我们这么个小的程序完全没必要那么做，那么，我们接下来就去 “白嫖” Deta Space的功能来完成后端的部署，这样我们的前端在任何地方都可以调取该API，前端就是个客户端APP而已。

### 我们先把Fastapi的业务逻辑完成

很简单，我们只需要在main.app中写就可以了：

**首先导入需要的模块：**
```python
from fastapi import FastAPI, File, UploadFile, Form, Depends
import pandas as pd
import os
```

因为我们读取表格文件，然后生成shell脚本这个任务的返回结果，是我们的一个依赖（一直要用的东西。

**所以我们先把这个依赖写出来，**
可以看到，依赖函数其实就是一个去掉了路由的路由函数！
```python
async def load_file(file: UploadFile = File(...), 
              sampleName: int = Form(...),
              barcode: int = Form(...),
              chip: int = Form(...),
              lane: int = Form(...),
              dataPath: int = Form(...))->pd.DataFrame:
    contents = await file.read()
    df = pd.read_excel(contents)
    df = df.fillna('')
    df = df.iloc[:,[sampleName-1,barcode-1,chip-1,lane-1,dataPath-1]]
    df.columns = ["sampleName", "barcode", "chip", "lane", "dataPath"]
    return df
```

**然后是第一个路由函数：**

就像前面说的，我们的冗余处理，简单的加后缀即可。
```python
##  将冗余列加后缀
def makenames(df:pd.DataFrame,col:str)->pd.DataFrame:
    s='_'+df.groupby(col).cumcount().add(1).astype(str)
    df.loc[:,col]+=s.mask(s=="_1","")
    return df
    
    
@app.post("/upload")
async def process(path:str = Form(...),df:pd.DataFrame = Depends(load_file)):
    isBlank = 0
    isDuplicate = 0
    remotePath = path
    ## 判断dataPath列书否为""
    if "" in df.dataPath.tolist():
        isBlank = 1
        df = df[df.dataPath!=""]
    if df.sampleName.duplicated().sum() > 0:
        isDuplicate = 1
        df = makenames(df,"sampleName")
    df = df.assign(filename1=df.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"1.fq.gz"]),axis=1))
    df = df.assign(filename2=df.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"2.fq.gz"]),axis=1))
    cmd1 = "\n".join(df.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename1']) + " " + os.path.join(remotePath,row['sampleName'] + "_R1.fastq.gz"),axis=1).tolist())
    cmd2 = "\n".join(df.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename2']) + " " + os.path.join(remotePath,row['sampleName'] + "_R2.fastq.gz"),axis=1).tolist())
    return {"rawData":df.to_dict(orient="records"),"cmd1":cmd1,"cmd2":cmd2,"HasBlank":isBlank,"hasDup":isDuplicate}
```


**最后是第二个函数**
通过提交样本名对应文件，修改样本名为真实数据。

```python
@app.post("/changenames")
async def change(file2:UploadFile= File(...),path:str = Form(...),df2:pd.DataFrame = Depends(load_file)):
    remotePath = path
    isDuplicate = 0
    df = pd.read_csv(file2.file,header=None,delimiter="\t")
    df.columns = ["sampleName","trueName"]
    df2 = df2.merge(df,how="left",on="sampleName")
    df2 = df2.drop(columns=["sampleName"])
    df2 = df2.dropna()
    df2 = df2.rename(columns={"trueName":"sampleName"})
    if df.sampleName.duplicated().sum() > 0:
        isDuplicate = 1
        df2 = makenames(df,"sampleName")
    df2 = df2.assign(filename1=df2.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"1.fq.gz"]),axis=1))
    df2 = df2.assign(filename2=df2.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"2.fq.gz"]),axis=1))
    cmd1 = "\n".join(df2.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename1']) + " " + os.path.join(remotePath,row['sampleName'] + "_R1.fastq.gz"),axis=1).tolist())
    cmd2 = "\n".join(df2.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename2']) + " " + os.path.join(remotePath,row['sampleName'] + "_R2.fastq.gz"),axis=1).tolist())
    return {"rawData":df2.to_dict(orient="records"),"cmd1":cmd1,"cmd2":cmd2,"HasBlank":0,"hasDup":isDuplicate}
```

**跨域CORS问题：**
别忘了，我们要随处访问该API，所以跨域问题需要暴力移除：
添加几行代码在头部：

```python
from fastapi.middleware.cors import CORSMiddleware

## 1. Create the FastAPI instance
app = FastAPI()

## 2. Add the CORS middleware
origins = [
    "*"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### 注册Deta Space帐户

创建免费的Deta Space帐户
接下来，创建一个免费帐户[Deta Space](https://deta.space/signup?dev_mode=true&ref=fastapi)，您只需要一个电子邮件和密码。

您甚至不需要信用卡，**但请确保在注册时启用开发者模式**。


### 安装space cli工具

```bash
curl -fsSL https://get.deta.dev/space-cli.sh | sh
```

### 登陆cli

为了使用 Deta Space 验证您的 CLI，您将需要一个访问令牌。

要获取此令牌，请打开您的Deta Space canvas，打开Teletype（画布底部的命令栏），然后单击“设置”。从那里，选择生成令牌并复制生成的令牌。

![token](https://picshack.net/ib/J5YNdUOajB.png)

现在从 Space CLI 运行`space login`。将令牌粘贴到 CLI 提示符并按 Enter 键后，您应该会看到一条确认消息。

```bash
$ space login

To authenticate the Space CLI with your Space account, generate a new access token in your Space settings and paste it below:

*****************************************

👍 Login Successful!
```

### 在空间中创建新的工程并提交

#### 创建工程
现在您已经使用 Space CLI 进行了身份验证，可以使用它来创建一个新的工程:

```bash
$ space new

# What is your project's name? >$ Glims2Excel

```
Space CLI 会要求您为项目命名，我们将称之为*Glims2Excel*。

然后，它会尝试自动检测您正在使用的框架或语言，并向您展示它找到的内容。在我们的例子中，它将通过以下消息识别 Python 应用程序，提示您确认：

```python
⚙️ No Spacefile found, trying to auto-detect configuration ...
👇 Deta detected the following configuration:

Micros:
name: Glims2Excel
 L src: .
 L engine: python3.9

# Do you want to bootstrap "Glims2Excel" with this configuration? (y/n)$ y

```

确认后，您的项目将在 Deta Space 中创建。[Builder](https://deta.space/builder) 是一个工具箱，可帮助您在 Deta Space 中创建和管理应用程序。

CLI 还将在本地Glims2Excel目录中创建一个Spacefile。Spacefile是一个配置文件，告诉 Deta Space 如何运行您的应用程序。您的应用程序的内容Spacefile如下：

```yaml
v: 0
micros:
  - name: Glims2Excel
    src: .
    engine: python3.9
```

它是一个yaml文件，您可以使用它来添加计划任务等功能或修改应用程序的运行方式，我们稍后将执行此操作。

#### 在Spacefile中定义运行命令

Spacefile 中的命令告诉 Space 应执行什么命令来启动您的应用程序。在这种情况下，它应该是uvicorn main:app。注意，我还添加了`public_routes`选项，目的是令我的API可以公开访问。

```yaml
v: 0
micros:
  - name: Glims2Excel
    src: .
    engine: python3.9
    run: uvicorn main:app
    public_routes:
      - "/*"
```

#### 添加requirements.txt

在正式上传之前，务必添加依赖文件，该依赖文件应该和main.py在同一个位置，同一个目录下，否则无法运行：

```
# requirements.txt
fastapi
uvicorn[standard]
python-multipart ##处理Formdata
pandas ##分析表格数据
openpyxl ## 读取excel
```

#### 正式部署

运行`space push`即可。你可以在你的Canvas中获得一个app，本文的app地址为 [https://glims2excel-1-w4936186.deta.app/docs](https://glims2excel-1-w4936186.deta.app/docs) 。

### 总结

经过以上步骤，我们成功的将自己的api部署到了容器中，我们可以在任何地方采用任何方式来和该app交互。下一篇文章，我们将把前端页面写出来。