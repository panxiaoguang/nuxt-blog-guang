<template>
  <article class="group border m-2 rounded-2xl overflow-hidden shadow-lg text-zinc-700">
    <NuxtLink to="/">
      <div
        class="lg:h-48 md:h-36 w-full object-cover object-center group-hover:scale-[1.05] transition-all duration-500"
      >
        <LogoConfused />
      </div>
      <div class="p-5">
        <h2 class="text-3xl font-semibold text-black pb-1 group-hover:text-sky-700">
          无可用文档
        </h2>

        <div class="flex group-hover:underline text-sky-700 items-center pt-2">
          <p>回到首页</p>
          <LogoArrow />
        </div>
      </div>
    </NuxtLink>
  </article>
</template>
