<script setup>
import { Share } from 'vue3-social-share';
import 'vue3-social-share/lib/index.css'
const baseUrl = 'https://www.xiaohanyuese.top';
defineProps(['path', 'title', 'description', 'image']);
</script>

<template>
    <Share :url="baseUrl + path" :title="title" source="modifySource" :description="description"
        :image="image" :google="true"></Share>
</template>