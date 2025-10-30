"""
Security Middleware for Kanad Platform
Implements security headers and protections
"""

from fastapi import Request
from fastapi.responses import Response
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.types import ASGIApp
import os


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """
    Adds security headers to all responses
    Implements OWASP recommended headers
    """

    def __init__(self, app: ASGIApp):
        super().__init__(app)
        self.environment = os.getenv("ENVIRONMENT", "development")

    async def dispatch(self, request: Request, call_next):
        response = await call_next(request)

        # Content Security Policy (CSP)
        # Prevents XSS attacks by controlling resource loading
        if self.environment == "production":
            csp = "; ".join([
                "default-src 'self'",
                "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://3dmol.org https://accounts.google.com https://www.gstatic.com",
                "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com",
                "font-src 'self' https://fonts.gstatic.com",
                "img-src 'self' data: https: blob:",
                "connect-src 'self' https://pubchem.ncbi.nlm.nih.gov https://cactus.nci.nih.gov https://accounts.google.com",
                "frame-src 'self' https://accounts.google.com",
                "object-src 'none'",
                "base-uri 'self'",
                "form-action 'self'",
                "frame-ancestors 'none'",
                "upgrade-insecure-requests",
            ])
        else:
            # More permissive CSP for development
            csp = "; ".join([
                "default-src 'self' 'unsafe-inline' 'unsafe-eval'",
                "connect-src 'self' http://localhost:* https://pubchem.ncbi.nlm.nih.gov https://cactus.nci.nih.gov",
                "img-src 'self' data: https: blob:",
            ])

        response.headers["Content-Security-Policy"] = csp

        # X-Content-Type-Options
        # Prevents MIME type sniffing
        response.headers["X-Content-Type-Options"] = "nosniff"

        # X-Frame-Options
        # Prevents clickjacking attacks
        response.headers["X-Frame-Options"] = "DENY"

        # X-XSS-Protection
        # Legacy XSS protection (modern browsers use CSP)
        response.headers["X-XSS-Protection"] = "1; mode=block"

        # Strict-Transport-Security (HSTS)
        # Forces HTTPS (only in production)
        if self.environment == "production":
            response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains; preload"

        # Referrer-Policy
        # Controls referrer information
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"

        # Permissions-Policy (formerly Feature-Policy)
        # Controls browser features
        permissions = "; ".join([
            "geolocation=()",
            "microphone=()",
            "camera=()",
            "payment=()",
            "usb=()",
            "magnetometer=()",
            "gyroscope=()",
            "accelerometer=()",
        ])
        response.headers["Permissions-Policy"] = permissions

        # X-Permitted-Cross-Domain-Policies
        # Prevents Adobe Flash and PDF cross-domain requests
        response.headers["X-Permitted-Cross-Domain-Policies"] = "none"

        # Cache-Control for sensitive endpoints
        if "/api/auth/" in request.url.path or "/api/admin/" in request.url.path:
            response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, private"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"

        return response


class RequestLoggingMiddleware(BaseHTTPMiddleware):
    """
    Logs all requests for security auditing
    """

    async def dispatch(self, request: Request, call_next):
        # Log request
        client_host = request.client.host if request.client else "unknown"
        print(f"→ {request.method} {request.url.path} from {client_host}")

        # Process request
        response = await call_next(request)

        # Log response status
        status_emoji = "✓" if response.status_code < 400 else "✗"
        print(f"{status_emoji} {request.method} {request.url.path} → {response.status_code}")

        return response
