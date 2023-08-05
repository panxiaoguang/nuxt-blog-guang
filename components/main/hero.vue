<script setup>
import aplayer from 'vue-aplayer-next'
</script>
<template>
  <div class="container mx-auto">
    <div class="grid grid-cols-1 sm:grid-cols-2 items-center">
      <div class="px-6">
        <h1 class="text-black font-semibold leading-tight text-4xl md:text-5xl my-5">
          欢迎访问我的博客
        </h1>
        <p>
          分享，生物信息学和编程经验，记录生活点滴
        </p>
        <ClientOnly>
          <div class="mt-5">
            <aplayer autoplay :music="{
              title: 'Wish',
              artist: ' - S.E.N.S.(神思者)',
              src: '/music/Wish.mp3',
              pic: 'https://picshack.net/ib/7HhaUKjMg1.jpg'
            }" />
          </div>
        </ClientOnly>
    </div>
    <div class="px-6 justify-self-center">
      <LogoDog />
    </div>
  </div>
</div></template>
