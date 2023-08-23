<script setup>
useHead({
  title: '主页',
  meta: [
    {
      name: 'description',
      content:
        '潘晓光的个人博客，分享，生物信息学和编程经验，记录生活点滴',
    },
    {
      name: 'keywords',
      content:
        '生信, 生物信息, 博客, 基因组, 转录组, 宏基因组, 编程, pytho, R, 生物信息学',
    },
    {
      name: 'baidu-site-verification',
      content:
        'codeva-ZVa7RA3e6R',
    },
  ],
  titleTemplate: '晓寒月色 - %s',
})

</script>

<template>
  <main class="container max-w-5xl mx-auto text-zinc-600">
    <MainHero />
    <MainRecent />
    <MainTrending />
  </main>
</template>
