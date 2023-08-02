<script setup lang="ts">
import type { BlogPost } from '@/types/blog'

const { path } = useRoute()
const articles = await queryContent(path).findOne()
const tocs = articles.body.toc
const data = computed<BlogPost>(() => {
  return {
    title: articles.title || 'no-title available',
    description: articles.description || 'no-descriptoin available',
    image: articles.image || '/nuxt-blog/no-image_cyyits.png',
    alt: articles.alt || 'no alter data available',
    ogImage: articles.ogImage || '/nuxt-blog/no-image_cyyits.png',
    date: articles.date || 'not-date-available',
    tags: articles.tags || [],
    published: articles.published || false,
  }
})

useHead({
  title: data.value.title || '',
  meta: [
    { name: 'description', content: data.value.description },
    {
      name: 'description',
      content: data.value.description,
    },
    // Test on: https://developers.facebook.com/tools/debug/ or https://socialsharepreview.com/
    { property: 'og:site_name', content: '晓寒月色' },
    { hid: 'og:type', property: 'og:type', content: 'website' },
    {
      property: 'og:title',
      content: data.value.title,
    },
    {
      property: 'og:description',
      content: data.value.description,
    },
    {
      property: 'og:image',
      content: data.value.ogImage || data.value.image,
    },
  ],
})
</script>

<template>
  <main class="px-6 container mx-auto">
    <header class="max-w-5xl mx-auto">
      <h1 class="text-xl md:text-3xl lg:text-4xl m-7 font-bold text-center">
        {{ data.title || '' }}
      </h1>
      <NuxtImg :src="data.image || ''" :alt="data.alt || ''"
        class="m-auto rounded-2xl shadow-lg h-52 md:h-96 w-4/5 content-center object-cover" />
      <p class="text-xs sm:text-sm my-3 max-w-3xl mx-auto text-center text-zinc-600">
        {{ data.description }}
      </p>
      <div class="flex w-full justify-center text-xs md:text-base my-8">
        <div class="md:flex text-black content-center gap-8 text-xs sm:text-sm">
          <div class="flex items-center">
            <LogoDate />
            <p>{{ data.date || '' }}</p>
          </div>
          <div class="flex items-center gap-2 flex-wrap my-5">
            <LogoTag />
            <template v-for="tag in data.tags" :key="tag">
              <span class="bg-gray-200 rounded-md px-2 py-1">{{ tag }}</span>
            </template>
          </div>
        </div>
      </div>
    </header>
    <div class="sm:grid grid-cols-12 max-w-7xl container">
      <div
        class="col-span-12 prose prose-pre:max-w-xs sm:prose-pre:max-w-full prose-sm sm:prose-base md:prose-lg prose-h1:no-underline  mx-auto prose-zinc prose-img:rounded-lg">
        <ContentRenderer :value="articles">
          <template #empty>
            <p>No content found.</p>
          </template>
        </ContentRenderer>
      </div>
      <div class="fixed top-36 left-16 max-w-[200px] max-h-80 overflow-y-scroll hidden sm:block">
        <div class="pb-2" v-if="tocs && tocs.links">
          <h3 class="text-lg">文章索引：</h3>
        </div>
        <ul v-if="tocs && tocs.links" class="list-disc">
          <li v-for="link in tocs.links" :key="link.text">
            <a :href="`#${link.id}`" class="text-slate-400 hover:text-sky-400">
              {{ link.text }}
            </a>
          </li>
        </ul>
      </div>
    </div>
  </main>
</template>

