<script setup lang="ts">

</script>

<template>
  <div class="py-5">
    <div class="sm:grid grid-cols-10 px-6 py-5 sm:py-9 gap-5 container max-w-5xl  mx-auto">
      <div class="col-span-5  max-w-md">
        <div class="flex justify-between">
          <div>
            <h1 class="text-xl sm:text-4xl  pb-2 font-bold ">
              小光amateur
            </h1>
            <div class="my-3 space-x-3 pb-10">
              <NuxtLink to="https://github.com/panxiaoguang" target="_blank" class="p-2 bg-gray-300 rounded-md"
                aria-label="Github">
                <Icon name="fa:github" size="1em" />
              </NuxtLink>
              <NuxtLink to="https://space.bilibili.com/6078012/" target="_blank" class="p-2 bg-gray-300 rounded-md"
                aria-label="BiliBili">
                <Icon name="fa6-brands:bilibili" size="1em" />
              </NuxtLink>
              <NuxtLink to="https://www.jianshu.com/u/92d23bce2ef3" target="_blank" class="p-2 bg-gray-300 rounded-md"
                aria-label="Jianshu">
                <Icon name="mdi:emoticon-happy" size="1em" />
              </NuxtLink>
              <NuxtLink to="https://www.zhihu.com/people/luo-tian-bao-92" target="_blank"
                class="p-2 bg-gray-300 rounded-md" aria-label="Zhihu">
                <Icon name="tabler:brand-zhihu" size="1em" />
              </NuxtLink>
            </div>
          </div>
          <div class="sm:hidden block col-span-5 pb-5">
            <NuxtImg src="/avatar.jpg" width="125" height="115" quality="50" class="rounded-md" />
          </div>
        </div>
        <h3 class="text-base sm:text-3xl pb-7 sm:pb-12">
          哥本哈根生物学博士，华大基因研究员生物信息副研。
        </h3>

        <p>
          潘晓光，第一批90后，已在社会上摸爬滚打了4年有余……</p>
        <p>本硕就读于山东省内的“985工程”，“211工程”的山东大学，毕业后不辍学业，于2018年7月进入青岛华大基因研究院，继续在生物专业领域深耕科研。主攻基因组编辑生物信息学分析，次要方向为染色体外环状DNA
          (eccDNA) 在癌症等多种疾病中的潜在作用及标记基因发现。</p>
        <p>2022年11月1日，在公司联合培养项目下继续进入丹麦哥本哈根大学深造，攻读博士研究生学位，方向为基因组编辑。</p>
      </div>
      <div class="hidden sm:block col-span-5">
        <h3 class="text-base sm:pb-12 sm:text-3xl">
          发表文献
        </h3>
        <ol class="list-decimal list-inside">
          <li>Pan, X., Qu, K., Yuan, H., Xiang, X., Anthon, C., Pashkova, L., Liang, X., Han, P., Corsi, G. I., & Xu, F. (2022). Massively targeted evaluation of therapeutic CRISPR off-targets in cells. Nature communications, 13(1), 4049.</li>
          <li>Corsi, G. I., Qu, K., Alkan, F., Pan, X., Luo, Y., & Gorodkin, J. (2022). CRISPR/Cas9 gRNA activity depends on free energy changes and on the target PAM context. Nature communications, 13(1), 3006.</li>
          <li>Jiang, X., Pan, X., Li, W., Han, P., Yu, J., Li, J., Zhang, H., Lv, W., Zhang, Y., & He, Y. (2023). Genome-wide characterization of extrachromosomal circular DNA in gastric cancer and its potential role in carcinogenesis and cancer progression.</li>
          <li>Lv, W., Pan, X., Han, P., Wang, Z., Feng, W., Xing, X., Wang, Q., Qu, K., Zeng, Y., & Zhang, C. (2022). Circle‐Seq reveals genomic and disease‐specific hallmarks in urinary cell‐free extrachromosomal circular DNAs. Clinical and Translational Medicine, 12(4), e817.</li>
        </ol>
      </div>
    </div>
  </div>
</template>
