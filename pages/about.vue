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
        <h3 class="text-base sm:text-3xl font-semibold pb-7 sm:pb-12">
          哥本哈根生物学博士，华大基因研究员生物信息副研。
        </h3>

        <p>
          潘晓光，第一批90后，已在社会上摸爬滚打了4年有余……</p>
        <p>本硕就读于山东省内的“985工程”，“211工程”的<span class="text-orange-600">山东大学</span>，毕业后不辍学业，于2018年7月进入<span class="text-orange-600">青岛华大基因研究院</span>，继续在生物专业领域深耕科研。主攻<span class="text-orange-600">基因组编辑生物信息学分析，次要方向为染色体外环状DNA</span>
          (eccDNA) 在癌症等多种疾病中的潜在作用及标记基因发现。</p>
        <p>2022年11月1日，在公司联合培养项目下继续进入<span class="text-orange-600">丹麦哥本哈根大学</span>深造，攻读博士研究生学位，方向为基因组编辑。</p>
      </div>
      <div class="hidden sm:block col-span-5">
        <h3 class="text-base sm:pb-12 sm:text-3xl">
          发表文献
        </h3>
        <ol class="list-decimal list-inside">
          <li><a class="text-sky-400 after:content-['_↗']" href="https://orcid.org/0000-0002-2688-9245" target="_blank">Pan, X.</a>, Qu, K., Yuan, H., Xiang, X., Anthon, C., Pashkova, L., Liang, X., Han, P., Corsi, G. I., & Xu, F. (2022). Massively targeted evaluation of therapeutic CRISPR off-targets in cells. <span class="italic text-orange-600">Nature communications</span>, 13(1), 4049.</li>
          <li>Corsi, G. I., Qu, K., Alkan, F., <a class="text-sky-400 after:content-['_↗']" href="https://orcid.org/0000-0002-2688-9245" target="_blank">Pan, X.</a>, Luo, Y., & Gorodkin, J. (2022). CRISPR/Cas9 gRNA activity depends on free energy changes and on the target PAM context. <span class="italic text-orange-600">Nature communications</span>, 13(1), 3006.</li>
          <li>Jiang, X., <a class="text-sky-400 after:content-['_↗']" href="https://orcid.org/0000-0002-2688-9245" target="_blank">Pan, X.</a>, Li, W. et al. Genome-wide characterization of extrachromosomal circular DNA in gastric cancer and its potential role in carcinogenesis and cancer progression. <span class="italic text-orange-600">Cell. Mol. Life Sci.</span> 80, 191 (2023). </li>
          <li>Lv, W., <a class="text-sky-400 after:content-['_↗']" href="https://orcid.org/0000-0002-2688-9245" target="_blank">Pan, X.</a>, Han, P., Wang, Z., Feng, W., Xing, X., Wang, Q., Qu, K., Zeng, Y., & Zhang, C. (2022). Circle-Seq reveals genomic and disease-specific hallmarks in urinary cell-free extrachromosomal circular DNAs. <span class="italic text-orange-600">Clinical and Translational Medicine</span>, 12(4), e817.</li>
        </ol>
      </div>
    </div>
  </div>
</template>
