// https://nuxt.com/docs/api/configuration/nuxt-config
export default defineNuxtConfig({
  app: {
    head: {
      charset: 'utf-16',
      viewport: 'width=device-width,initial-scale=1',
      title: '晓寒月色',
      titleTemplate: '%s - 晓寒月色',
      meta: [{ name: 'description', content: '生物信息技术博客' }],
    },
    pageTransition: { name: 'page', mode: 'out-in' },
    layoutTransition: { name: 'layout', mode: 'out-in' },
  },
  sitemap: {
    strictNuxtContentPaths: true,
  },
  site: {
    url: 'https://www.xiaohanyuese.top',
  },

  typescript: {
    strict: true,
  },

  nitro: {
    prerender: {
      crawlLinks: true,
      routes: [
        '/',
      ],
    },
  },

  modules: [
    '@nuxt/content',
    '@nuxtjs/tailwindcss',
    'nuxt-icon',
    '@nuxt/image-edge',
    '@nuxtjs/robots',
    '@nuxtjs/fontaine',
    'nuxt-simple-sitemap',
  ],

  content: {
    highlight: {
      theme: 'dracula',
    },
    markdown: {
      remarkPlugins: [
        'remark-math',
      ],
      rehypePlugins: [
        'rehype-mathjax',
      ]
    }
  }
})
