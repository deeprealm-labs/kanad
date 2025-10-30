"""
Compute Resource Limits Middleware
Prevents abuse of computational resources by enforcing quotas and limits
"""
from fastapi import Request, HTTPException, status
from starlette.middleware.base import BaseHTTPMiddleware
from typing import Dict, Optional
from datetime import datetime, timedelta
from collections import defaultdict
import asyncio

# ============================================================================
# Configuration - Adjust these limits based on your Azure tier
# ============================================================================

# Maximum atoms per molecule (prevents huge molecules)
MAX_ATOMS_PER_MOLECULE = 50  # ~50 atoms is reasonable for VQE/SQD

# Maximum basis set size
ALLOWED_BASIS_SETS = ['sto-3g', '3-21g', '6-31g', '6-31g*', 'cc-pvdz']
MAX_BASIS_FUNCTIONS = 200  # Limits computational complexity

# Maximum iterations for optimization
MAX_VQE_ITERATIONS = 500
MAX_SQD_ITERATIONS = 100

# Maximum concurrent jobs per user
MAX_CONCURRENT_JOBS_PER_USER = 3

# Maximum experiments per user per hour (rate limiting)
MAX_EXPERIMENTS_PER_HOUR = 10

# Maximum experiments per user per day
MAX_EXPERIMENTS_PER_DAY = 50

# Backend restrictions
ALLOWED_BACKENDS_FREE_TIER = ['statevector', 'aer_simulator']
ALLOWED_BACKENDS_PREMIUM = ['statevector', 'aer_simulator', 'ibm_quantum', 'bluequbit']

# User tier system
class UserTier:
    FREE = "free"
    PREMIUM = "premium"
    ADMIN = "admin"

# ============================================================================
# In-Memory Rate Limiting (use Redis in production)
# ============================================================================

class RateLimiter:
    def __init__(self):
        self.user_experiments: Dict[int, list] = defaultdict(list)
        self.user_jobs: Dict[int, int] = defaultdict(int)
        self._lock = asyncio.Lock()

    async def check_hourly_limit(self, user_id: int, limit: int = MAX_EXPERIMENTS_PER_HOUR) -> bool:
        """Check if user has exceeded hourly experiment limit"""
        async with self._lock:
            now = datetime.utcnow()
            one_hour_ago = now - timedelta(hours=1)

            # Clean old entries
            self.user_experiments[user_id] = [
                ts for ts in self.user_experiments[user_id]
                if ts > one_hour_ago
            ]

            # Check limit
            if len(self.user_experiments[user_id]) >= limit:
                return False

            # Add new entry
            self.user_experiments[user_id].append(now)
            return True

    async def check_daily_limit(self, user_id: int, limit: int = MAX_EXPERIMENTS_PER_DAY) -> bool:
        """Check if user has exceeded daily experiment limit"""
        async with self._lock:
            now = datetime.utcnow()
            one_day_ago = now - timedelta(days=1)

            # Count experiments in last 24 hours
            recent_experiments = [
                ts for ts in self.user_experiments[user_id]
                if ts > one_day_ago
            ]

            return len(recent_experiments) < limit

    async def check_concurrent_jobs(self, user_id: int, limit: int = MAX_CONCURRENT_JOBS_PER_USER) -> bool:
        """Check if user has too many concurrent jobs"""
        async with self._lock:
            return self.user_jobs[user_id] < limit

    async def increment_jobs(self, user_id: int):
        """Increment concurrent job counter"""
        async with self._lock:
            self.user_jobs[user_id] += 1

    async def decrement_jobs(self, user_id: int):
        """Decrement concurrent job counter"""
        async with self._lock:
            if self.user_jobs[user_id] > 0:
                self.user_jobs[user_id] -= 1

# Global rate limiter instance
rate_limiter = RateLimiter()

# ============================================================================
# Validation Functions
# ============================================================================

def validate_molecule_size(molecule_data: dict) -> tuple[bool, Optional[str]]:
    """Validate molecule doesn't exceed size limits"""
    # Count atoms
    atom_count = 0

    if "atoms" in molecule_data and molecule_data["atoms"]:
        atom_count = len(molecule_data["atoms"])
    elif "smiles" in molecule_data:
        # Rough estimate from SMILES (each non-H heavy atom)
        smiles = molecule_data["smiles"]
        # Simple heuristic: count uppercase letters (roughly heavy atoms)
        atom_count = sum(1 for c in smiles if c.isupper())

    if atom_count > MAX_ATOMS_PER_MOLECULE:
        return False, f"Molecule too large: {atom_count} atoms (max: {MAX_ATOMS_PER_MOLECULE})"

    return True, None

def validate_basis_set(basis_set: str) -> tuple[bool, Optional[str]]:
    """Validate basis set is allowed"""
    if basis_set.lower() not in ALLOWED_BASIS_SETS:
        return False, f"Basis set '{basis_set}' not allowed. Allowed: {', '.join(ALLOWED_BASIS_SETS)}"

    return True, None

def validate_iterations(method: str, iterations: int) -> tuple[bool, Optional[str]]:
    """Validate iteration count is within limits"""
    if method.lower() in ['vqe', 'excited_states']:
        if iterations > MAX_VQE_ITERATIONS:
            return False, f"VQE iterations too high: {iterations} (max: {MAX_VQE_ITERATIONS})"
    elif method.lower() == 'sqd':
        if iterations > MAX_SQD_ITERATIONS:
            return False, f"SQD iterations too high: {iterations} (max: {MAX_SQD_ITERATIONS})"

    return True, None

def validate_backend(backend: str, user_tier: str = UserTier.FREE) -> tuple[bool, Optional[str]]:
    """Validate backend is allowed for user tier"""
    allowed = ALLOWED_BACKENDS_FREE_TIER if user_tier == UserTier.FREE else ALLOWED_BACKENDS_PREMIUM

    if backend.lower() not in allowed:
        return False, f"Backend '{backend}' not allowed for {user_tier} tier. Allowed: {', '.join(allowed)}"

    return True, None

# ============================================================================
# Compute Limits Middleware
# ============================================================================

class ComputeLimitsMiddleware(BaseHTTPMiddleware):
    """
    Middleware to enforce compute resource limits on experiment submissions
    """

    async def dispatch(self, request: Request, call_next):
        # Only check compute limits on experiment submission endpoints
        if request.url.path in ['/api/experiments/submit', '/api/campaigns/create']:
            # Get user from request state (set by auth middleware)
            user = getattr(request.state, 'user', None)

            if user:
                # Check rate limits
                if not await rate_limiter.check_hourly_limit(user.id):
                    raise HTTPException(
                        status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                        detail=f"Hourly experiment limit exceeded ({MAX_EXPERIMENTS_PER_HOUR} per hour). Please try again later."
                    )

                if not await rate_limiter.check_daily_limit(user.id):
                    raise HTTPException(
                        status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                        detail=f"Daily experiment limit exceeded ({MAX_EXPERIMENTS_PER_DAY} per day). Upgrade to premium for higher limits."
                    )

                if not await rate_limiter.check_concurrent_jobs(user.id):
                    raise HTTPException(
                        status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                        detail=f"Too many concurrent jobs ({MAX_CONCURRENT_JOBS_PER_USER} max). Wait for current jobs to complete."
                    )

        response = await call_next(request)
        return response

# ============================================================================
# Helper Functions for Route Validation
# ============================================================================

async def validate_experiment_config(config: dict, user_tier: str = UserTier.FREE) -> None:
    """
    Validate experiment configuration against compute limits
    Raises HTTPException if validation fails
    """
    errors = []

    # Validate molecule size
    if "molecule" in config:
        valid, error = validate_molecule_size(config["molecule"])
        if not valid:
            errors.append(error)

    # Validate basis set
    if "configuration" in config:
        conf = config["configuration"]

        # Basis set
        if "basis" in conf:
            valid, error = validate_basis_set(conf["basis"])
            if not valid:
                errors.append(error)

        # Backend
        if "backend" in conf:
            valid, error = validate_backend(conf["backend"], user_tier)
            if not valid:
                errors.append(error)

        # Iterations
        if "max_iterations" in conf:
            method = conf.get("method", "vqe")
            valid, error = validate_iterations(method, conf["max_iterations"])
            if not valid:
                errors.append(error)

    if errors:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail={
                "error": "Experiment configuration exceeds compute limits",
                "violations": errors,
                "limits": {
                    "max_atoms": MAX_ATOMS_PER_MOLECULE,
                    "allowed_basis_sets": ALLOWED_BASIS_SETS,
                    "max_vqe_iterations": MAX_VQE_ITERATIONS,
                    "max_sqd_iterations": MAX_SQD_ITERATIONS,
                    "allowed_backends": ALLOWED_BACKENDS_FREE_TIER if user_tier == UserTier.FREE else ALLOWED_BACKENDS_PREMIUM
                }
            }
        )

def get_user_tier(user) -> str:
    """Get user tier from user object"""
    if not user:
        return UserTier.FREE

    # Admin users have unlimited access
    if hasattr(user, 'role') and user.role.value == 'admin':
        return UserTier.ADMIN

    # Check if user has premium subscription (implement your logic)
    # For now, all non-admin users are free tier
    return UserTier.FREE

# ============================================================================
# Quota Information Endpoint
# ============================================================================

async def get_user_quota_info(user_id: int) -> dict:
    """Get user's current quota usage"""
    now = datetime.utcnow()
    one_hour_ago = now - timedelta(hours=1)
    one_day_ago = now - timedelta(days=1)

    # Get experiment history
    experiments = rate_limiter.user_experiments.get(user_id, [])

    hourly_count = len([ts for ts in experiments if ts > one_hour_ago])
    daily_count = len([ts for ts in experiments if ts > one_day_ago])
    concurrent_jobs = rate_limiter.user_jobs.get(user_id, 0)

    return {
        "experiments_this_hour": hourly_count,
        "hourly_limit": MAX_EXPERIMENTS_PER_HOUR,
        "experiments_today": daily_count,
        "daily_limit": MAX_EXPERIMENTS_PER_DAY,
        "concurrent_jobs": concurrent_jobs,
        "max_concurrent_jobs": MAX_CONCURRENT_JOBS_PER_USER,
        "compute_limits": {
            "max_atoms": MAX_ATOMS_PER_MOLECULE,
            "allowed_basis_sets": ALLOWED_BASIS_SETS,
            "max_vqe_iterations": MAX_VQE_ITERATIONS,
            "max_sqd_iterations": MAX_SQD_ITERATIONS
        }
    }
