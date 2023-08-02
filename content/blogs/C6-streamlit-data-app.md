---
title: 用streamlit搭建数据交互式app
date: 2ed Aug 2023
description: 大数据时代，数据交互越来越重要，而带有图形界面的web应用程序的需求也日益明显
image: /blogs-img/streamlit.webp
alt: 用streamlit搭建数据交互式app
ogImage: /blogs-img/streamlit.webp
tags: ['Python','前端']
published: true
---

### streamlit的有趣特点

- 所有的程序，只要是前端交互页面发生变动或者说交互，代码就会从头到尾执行一遍
- 提供了非常多数据交互的组件，每个组件都可以返回数值，用来和别的组件交流
- 有特殊的缓存系统，防止长时间运行的程序成为瓶颈
- 因为程序从头至尾的顺序执行，异步的支持较差

### 搭建一个小型demo

#### 主页1

![](https://picshack.net/ib/xDRJygkePJ.png)

#### 主页2

![](https://picshack.net/ib/vSEYy2f9WT.png)

#### 主要功能

下机数据的文件名字往往都是通过`barcode`区分的，而我们分析的时候需要按照样本为其命名，这样更方便，也不容易搞错。一个方便的方法便是用列表中的样本自动生成新的文件名，并自动提供数据拷贝指令。

这里有一个坑：如果在一个lane上测序量不足的情况下，会出现加测的情况，也就意味着，一个样本名会对应多个文件。换句话说，样本名列有冗余，需要我们注意。

所以搭建这个app的目的是：通过upload一个包含测序下机路径数据的excel表格，设置需要获取的信息所在的列索引，本app自动把提取后的数据展示在网页表格中。然后，在第二个tab中自动把数据所需要的拷贝指令生成出来。

这里对待冗余样本用了很简单的逻辑，既自动生成_1,_2...后缀。

#### 直接上代码

```python
### streamlit_app/main.py
import streamlit as st
import pandas as pd
import os 

st.set_page_config(
    page_title="Excel toolkits",
    layout="wide"
)

@st.cache_data
def read_data(upload_file):
    # Read the uploaded file into a dataframe
    df = pd.read_excel(upload_file)

    # Rename the columns
    df.columns = ["n"+str(i+1) for i in range(df.shape[1])]

    # Fill missing values with empty string
    df = df.fillna('')

    return df

def read_data2(upload_file):
    # Read the file into a Pandas dataframe
    duiying = pd.read_excel(upload_file)

    # Rename the columns
    duiying.columns = ["sampleName","trueName"]

    return duiying

# This code makes sure that all of the names in the dataframe are unique by adding a number at the end of each name that is not unique.
# The code also adds a number at the end of each name that is unique, but not the first one.

def makenames(df,col):
    s='_'+df.groupby(col).cumcount().add(1).astype(str)
    df.loc[:,col]+=s.mask(s=="_1","")
    return df


## Add a title
st.title("Excel tools for Glims")
## Split the page into 2 columns
col1, col2 = st.columns([0.35,0.65],gap="medium")

## 这个是左边的栏目
with col1:
    ###组件数据交互
    upload_data = st.file_uploader(label="Excel file from Glims...",
                     type=['xlsx','xls'])
    SM_index = st.number_input(label="Index for sampleName:",
                    min_value=1,
                    max_value=62,
                    value=4)
    BC_index = st.number_input(label="Index for barcode:",
                    min_value=1,
                    max_value=62,
                    value=23)
    CP_index = st.number_input(label="Index for chip:",
                    min_value=1,
                    max_value=62,
                    value=27)
    LN_index = st.number_input(label="Index for lane:",
                    min_value=1,
                    max_value=62,
                    value=28)
    Path_index = st.number_input(label="Index for DataPath:",
                    min_value=1,
                    max_value=62,
                    value=48)
    Remote_path = st.text_input(label="Remote path:",
                                value="10.2.100.1:/pakpox/pob8d1/",
                                max_chars=1000)
    isSubmit = st.button("Submit!")
    st.divider()
    placeholder = st.empty()
    upload_data2 = st.file_uploader(label="Input Barcode with sampleName:",
                     type=['xlsx','xls'])

### 这个是右边的栏目
with col2:
    ### 分两个tab
    tab1, tab2 = st.tabs(["ViewData", "Scipts"])
    tab21, tab22 = tab2.tabs(["Script1", "Script2"])
    ### container是容器占位符，作用是可以把后面运行的命令得到的结果放到一个容器里面渲染在前面去
    upcontent = tab1.container()
    upcontent.subheader("View of the data:")
    if upload_data is not None: ###这个判断很重要，不然会报错

        ###这是一个缓存函数
        dataframes = read_data(upload_data)
        df = dataframes.iloc[:,[SM_index-1,BC_index-1,CP_index-1,LN_index-1,Path_index-1]]
        df.columns = ["sampleName", "barcode", "chip", "lane", "dataPath"]
        ## if sampleName is duplicated, fix the name
        if df.sampleName.duplicated().any():
            tab1.warning("Duplicated sampleName found, fixing...")
        df = makenames(df,"sampleName")
        ### upload2的作用是，有的excel文件并没有样本名列，它提供了额外的样本名和barcode对应关系，需要我们提交后再次上传
        if upload_data2 is not None:
            placeholder.info("After upload data, please click the button again!")
            df2 = read_data2(upload_data2)
            df = df.merge(df2,how="left",on="sampleName")
            df = df.drop(columns=["sampleName"])
            df = df.rename(columns={"trueName":"sampleName"})
            df = df[~df["sampleName"].isna()]
        ### 这个判断导致延迟执行，每次更新数据后，只有重新点击提交按钮才会重新执行
        if isSubmit:
            upcontent.dataframe(df,hide_index=True,height=600)
            df = df.assign(filename1=df.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"1.fq.gz"]),axis=1))
            df = df.assign(filename2=df.apply(lambda row : "_".join([row['chip'],row['lane'],str(row['barcode']),"2.fq.gz"]),axis=1))
            cmd1 = df.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename1']) + " " + os.path.join(Remote_path,row['sampleName'] + "_R1.fastq.gz"),axis=1).tolist()
            cmd2 = df.apply(lambda row: "scp " + os.path.join(row['dataPath'],row['filename2']) + " " + os.path.join(Remote_path,row['sampleName'] + "_R2.fastq.gz"),axis=1).tolist()
            tab21.code("\n".join(cmd1))
            tab22.code("\n".join(cmd2))
```


### 总结

这个小demo用到了不少知识点

1. 函数前面加`@st.cache_data`可以起到缓存作用
2. 因为streamlit是从头到尾的执行，所以`if`判断很重要，
例如upload2其实从逻辑上是最后执行的，但是因为这个特性加在了代码前面
3. 通常我们加载数据后，都会有信息的提示，例如“数据加载成功”，“检测到冗余！”等。然后streamlit难以实现异步，导致要么该提示一闪而过，要么就永久留在页面上，不然会对后面运行的程序带来影响。所以
`container`可以当作占位符，把想要渲染的组件先于逻辑放在前端

### 下一篇预告

streamlit是一个快速搭建网页小程序的工具，但是我们看到，目前一个很简单的需求，实现起来都有重重限制。而且其异步的难以实现更是让用户体验非常糟糕。

下一篇，我们将使用基于**next.js**和**Chakra UI**的[reflex](https://reflex.dev/)来再次搭建该网页小程序。看看性能会有提升么？