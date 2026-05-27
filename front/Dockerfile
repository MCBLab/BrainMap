FROM node:22-alpine

WORKDIR /app

# Install dependencies first for caching
COPY package*.json ./
RUN npm install

# Copy the rest of the frontend code
COPY . .

# Expose Vite's default port
EXPOSE 5173

# Start the dev server and bind to all interfaces
CMD ["npm", "run", "dev", "--", "--host", "0.0.0.0"]
