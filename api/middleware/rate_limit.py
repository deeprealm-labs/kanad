"""
Rate Limiting Middleware for Kanad Platform
Prevents brute force attacks and API abuse
"""

from datetime import datetime, timedelta
from typing import Dict, Optional
from fastapi import Request, HTTPException, status
from collections import defaultdict
import asyncio

# In-memory rate limit storage (use Redis for production multi-server setup)
_rate_limit_storage: Dict[str, list] = defaultdict(list)
_cleanup_task: Optional[asyncio.Task] = None


class RateLimiter:
    """Rate limiter with sliding window algorithm"""

    def __init__(self, requests: int, window_seconds: int):
        """
        Initialize rate limiter

        Args:
            requests: Maximum number of requests allowed
            window_seconds: Time window in seconds
        """
        self.max_requests = requests
        self.window_seconds = window_seconds

    def _get_key(self, request: Request, identifier: str = None) -> str:
        """Get rate limit key for request"""
        if identifier:
            return f"custom:{identifier}"

        # Use IP address + endpoint as key
        ip = request.client.host if request.client else "unknown"
        endpoint = request.url.path
        return f"{ip}:{endpoint}"

    def _cleanup_old_requests(self, key: str, cutoff_time: datetime):
        """Remove requests older than window"""
        if key in _rate_limit_storage:
            _rate_limit_storage[key] = [
                ts for ts in _rate_limit_storage[key]
                if ts > cutoff_time
            ]

    async def check(self, request: Request, identifier: str = None) -> bool:
        """
        Check if request is allowed under rate limit

        Args:
            request: FastAPI request object
            identifier: Optional custom identifier (e.g., email)

        Returns:
            True if allowed, False if rate limited
        """
        key = self._get_key(request, identifier)
        now = datetime.utcnow()
        cutoff = now - timedelta(seconds=self.window_seconds)

        # Clean old requests
        self._cleanup_old_requests(key, cutoff)

        # Count recent requests
        recent_requests = len(_rate_limit_storage[key])

        if recent_requests >= self.max_requests:
            return False

        # Add current request
        _rate_limit_storage[key].append(now)
        return True

    async def get_remaining(self, request: Request, identifier: str = None) -> int:
        """Get remaining requests in current window"""
        key = self._get_key(request, identifier)
        now = datetime.utcnow()
        cutoff = now - timedelta(seconds=self.window_seconds)

        self._cleanup_old_requests(key, cutoff)
        recent_requests = len(_rate_limit_storage[key])

        return max(0, self.max_requests - recent_requests)

    async def get_reset_time(self, request: Request, identifier: str = None) -> Optional[datetime]:
        """Get time when rate limit resets"""
        key = self._get_key(request, identifier)

        if key not in _rate_limit_storage or not _rate_limit_storage[key]:
            return None

        oldest = min(_rate_limit_storage[key])
        return oldest + timedelta(seconds=self.window_seconds)


# Predefined rate limiters
login_limiter = RateLimiter(requests=5, window_seconds=900)  # 5 per 15 minutes
register_limiter = RateLimiter(requests=3, window_seconds=3600)  # 3 per hour
otp_limiter = RateLimiter(requests=3, window_seconds=300)  # 3 per 5 minutes
admin_limiter = RateLimiter(requests=100, window_seconds=60)  # 100 per minute


async def rate_limit_dependency(
    request: Request,
    limiter: RateLimiter,
    identifier: str = None,
    error_message: str = None
):
    """
    FastAPI dependency for rate limiting

    Args:
        request: FastAPI request
        limiter: RateLimiter instance to use
        identifier: Optional custom identifier
        error_message: Custom error message

    Raises:
        HTTPException: 429 if rate limited
    """
    allowed = await limiter.check(request, identifier)

    if not allowed:
        remaining = await limiter.get_remaining(request, identifier)
        reset_time = await limiter.get_reset_time(request, identifier)

        headers = {
            "X-RateLimit-Limit": str(limiter.max_requests),
            "X-RateLimit-Remaining": str(remaining),
        }

        if reset_time:
            headers["X-RateLimit-Reset"] = reset_time.isoformat()

        message = error_message or f"Rate limit exceeded. Try again in {limiter.window_seconds} seconds."

        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=message,
            headers=headers
        )


# Background cleanup task
async def cleanup_rate_limits():
    """Periodically clean up old rate limit entries"""
    while True:
        await asyncio.sleep(300)  # Run every 5 minutes

        now = datetime.utcnow()
        keys_to_delete = []

        for key, timestamps in _rate_limit_storage.items():
            # Remove timestamps older than 1 hour
            cutoff = now - timedelta(hours=1)
            _rate_limit_storage[key] = [ts for ts in timestamps if ts > cutoff]

            # Mark empty keys for deletion
            if not _rate_limit_storage[key]:
                keys_to_delete.append(key)

        # Delete empty keys
        for key in keys_to_delete:
            del _rate_limit_storage[key]

        if keys_to_delete:
            print(f"✓ Rate limit cleanup: removed {len(keys_to_delete)} expired entries")


def start_cleanup_task():
    """Start background cleanup task"""
    global _cleanup_task
    if _cleanup_task is None or _cleanup_task.done():
        _cleanup_task = asyncio.create_task(cleanup_rate_limits())
        print("✓ Rate limit cleanup task started")


def stop_cleanup_task():
    """Stop background cleanup task"""
    global _cleanup_task
    if _cleanup_task and not _cleanup_task.done():
        _cleanup_task.cancel()
        print("✓ Rate limit cleanup task stopped")
