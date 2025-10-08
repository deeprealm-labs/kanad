# Multi-stage Dockerfile for Kanad Backend

FROM python:3.11-slim as base

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libopenblas-dev \
    liblapack-dev \
    gfortran \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Install Kanad framework (assuming it's in parent directory)
# In production, this would be installed from PyPI or a separate package
RUN pip install -e /app/../kanad

# Create non-root user
RUN useradd -m -u 1000 kanad && chown -R kanad:kanad /app
USER kanad

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

# Run FastAPI with Uvicorn
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
