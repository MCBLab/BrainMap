import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  // O site e publicado em https://mcblab.github.io/BrainMap/, entao os assets
  // precisam ser referenciados a partir de /BrainMap/ e nao da raiz.
  base: '/BrainMap/',
})
