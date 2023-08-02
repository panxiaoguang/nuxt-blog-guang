<template>
  <div>
    <p class="text-black text-2xl py-1">
      小光amateur
    </p>
    <p class="py-2">
      一个生物信息学爱好者，喜欢编程，喜欢分享，喜欢记录生活点滴
    </p>
    <div class="my-3 space-x-3 pb-3">
      <NuxtLink
        to="https://github.com/panxiaoguang"
        target="_blank"
        class="p-2 bg-gray-300 rounded-md"
        aria-label="Github"
      >
        <Icon name="fa:github" size="1em" />
      </NuxtLink>
      <NuxtLink
        to="https://space.bilibili.com/6078012/"
        target="_blank"
        class="p-2 bg-gray-300 rounded-md"
        aria-label="BiliBili"
      >
        <Icon name="fa6-brands:bilibili" size="1em" />
      </NuxtLink>
      <NuxtLink
        to="https://www.jianshu.com/u/92d23bce2ef3"
        target="_blank"
        class="p-2 bg-gray-300 rounded-md"
        aria-label="Jianshu"
      >
        <Icon name="mdi:emoticon-happy" size="1em" />
      </NuxtLink>
      <NuxtLink
        to="https://www.zhihu.com/people/luo-tian-bao-92"
        target="_blank"
        class="p-2 bg-gray-300 rounded-md"
        aria-label="Zhihu"
      >
        <Icon name="tabler:brand-zhihu" size="1em" />
      </NuxtLink>
    </div>
  </div>
</template>
