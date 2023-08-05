// 引入 vssue
import Vssue from '@liamrad/vssue-vue3'
// 引入对应平台的 api 包
import GithubV4 from '@vssue/api-github-v4'


export default defineNuxtPlugin((nuxtApp) => {
    const runtimeConfig = useRuntimeConfig()
    // Doing something with nuxtApp
    nuxtApp.vueApp.use(Vssue, {
        // 设置要使用的平台 api
        api: GithubV4,

        // 在这里设置你使用的平台的 OAuth App 配置
        owner: 'panxiaoguang',
        repo: 'panxiaoguang.github.io',
        clientId: runtimeConfig.public.clientId,
        clientSecret: runtimeConfig.public.clientSecret, // 只有在使用某些平台时需要
    })
})