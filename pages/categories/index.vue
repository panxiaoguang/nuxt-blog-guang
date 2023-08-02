<script lang="ts" setup>
import { makeFirstCharUpper } from '@/utils/helper'

const topics = [
  '生物信息',
  '基因编辑',
  '单细胞测序',
  '转录组',
  '基因组',
  'R语言',
  'Python',
  'Julia',
  '生信软件',
  '随笔',
  '文献阅读',
  '前端',
  '后端',
  '机器学习',
]

useHead({
  title: '分类',
  meta: [
    {
      name: 'description',
      content:
        '这里是所有的主题分类',
    },
  ],
  titleTemplate: '晓寒月色 - %s',
})
</script>

<template>
  <main class="container max-w-5xl mx-auto text-zinc-600">
    <CategoryHero />
    <div class="flex flex-wrap px-6 mt-12 gap-3">
      <template v-for="topic in topics" :key="topic">
        <CategoryCard :title="makeFirstCharUpper(topic)" />
      </template>
    </div>
  </main>
</template>
