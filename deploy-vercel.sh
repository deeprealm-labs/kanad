#!/bin/bash
# Kanad Frontend - Vercel Deployment Script

set -e

echo "========================================="
echo "  Kanad Frontend Deployment to Vercel"
echo "========================================="

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

AZURE_API_URL="http://172.171.222.16/api"

echo ""
echo "${BLUE}Step 1: Checking prerequisites...${NC}"

# Check if we're in the right directory
if [ ! -d "web" ]; then
    echo "❌ Error: web directory not found"
    echo "Please run this script from the kanad root directory"
    exit 1
fi

# Check if vercel is installed
if ! command -v vercel &> /dev/null; then
    echo "${YELLOW}⚠️  Vercel CLI not installed${NC}"
    echo "Installing Vercel CLI..."
    npm install -g vercel
fi

# Check if logged into Vercel
if ! vercel whoami &> /dev/null; then
    echo "${YELLOW}⚠️  Not logged into Vercel${NC}"
    echo "Please login to Vercel:"
    vercel login
fi

echo "${GREEN}✓ Prerequisites OK${NC}"

echo ""
echo "${BLUE}Step 2: Preparing frontend...${NC}"

cd web

# Install dependencies
echo "Installing dependencies..."
npm install

# Create production environment file
echo "Creating production environment..."
cat > .env.production << EOF
# Kanad Web Application - Production Environment
NEXT_PUBLIC_API_URL=${AZURE_API_URL}
EOF

echo "${GREEN}✓ Frontend prepared${NC}"

echo ""
echo "${BLUE}Step 3: Testing build locally...${NC}"

# Build the project
npm run build

echo "${GREEN}✓ Build successful${NC}"

echo ""
echo "${BLUE}Step 4: Deploying to Vercel...${NC}"

# Deploy to Vercel
echo "Deploying to Vercel..."
vercel --prod --env NEXT_PUBLIC_API_URL="${AZURE_API_URL}"

echo ""
echo "${GREEN}=========================================${NC}"
echo "${GREEN}  Deployment Complete! 🚀${NC}"
echo "${GREEN}=========================================${NC}"
echo ""
echo "Next steps:"
echo "1. Update CORS on Azure backend:"
echo "   ${YELLOW}ssh kanadmin@172.171.222.16${NC}"
echo "   ${YELLOW}cd /opt/kanad && nano .env${NC}"
echo "   Add your Vercel URL to CORS_ORIGINS"
echo ""
echo "2. Test your deployment:"
echo "   - Visit your Vercel URL"
echo "   - Test API connectivity"
echo "   - Run quantum calculations"
echo ""
echo "3. (Optional) Set up custom domain on Vercel"
echo ""

cd ..

echo "Deployment URL will be shown above ⬆️"
echo ""
