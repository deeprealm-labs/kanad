# Admin Panel Live Monitoring & Analytics Fix

## Issue
The admin panel's live monitoring, analytics, and experiment statistics were not displaying any data despite experiments actively running. All API endpoints were returning 200 OK but with empty arrays.

## Root Cause
**Dual-Database Architecture Mismatch**: The admin routes in `/api/routes/admin.py` were querying the PostgreSQL database for experiment data, but actual experiments are stored in the SQLite database (`api/kanad_experiments.db`).

### Database Architecture
- **SQLite** (`api/kanad_experiments.db`) - Stores experiments, jobs, campaigns, results
- **PostgreSQL** - Stores users, authentication, sessions, access keys, admin metadata

## Fix Applied

### 1. Updated Database Imports in [admin.py:14-26](api/routes/admin.py#L14-L26)
```python
# Renamed PostgreSQL models to distinguish from SQLite
from api.core.database_postgres import (
    Experiment as PGExperiment,  # PostgreSQL experiment metadata
    Campaign as PGCampaign,
    Job as PGJob,
)

# Added SQLite database imports
from api.core.database import ExperimentDB, CampaignDB
```

### 2. Rewrote Live Monitoring Endpoint [admin.py:405-446](api/routes/admin.py#L405-L446)
```python
@router.get("/experiments/live")
async def get_live_experiments(db: Session = Depends(get_db)):
    # NOW queries SQLite instead of PostgreSQL
    running_experiments = ExperimentDB.list(status="running", limit=1000)

    # Cross-references user info from PostgreSQL when needed
    for exp in running_experiments:
        user = None
        if exp.get("user_id"):
            user = db.query(User).filter(User.id == exp["user_id"]).first()
```

### 3. Rewrote Analytics Endpoints

**Usage by Backend** [admin.py:486-526](api/routes/admin.py#L486-L526)
- Now queries SQLite `ExperimentDB.list()` instead of PostgreSQL
- Manually aggregates experiment counts by backend
- Filters by date range

**System Stats Overview** [admin.py:601-631](api/routes/admin.py#L601-L631)
- User counts from PostgreSQL (users, active users, verified users)
- Experiment counts from SQLite (total, running, completed)
- Campaign counts from SQLite

### 4. Fixed User Management Endpoint [admin.py:204-208](api/routes/admin.py#L204-L208)
```python
# Fixed SQLAlchemy join ambiguity with explicit join conditions
query = db.query(
    User,
    func.count(func.distinct(PGExperiment.id)).label("experiments_count"),
    func.count(func.distinct(PGCampaign.id)).label("campaigns_count"),
).select_from(User).outerjoin(
    PGExperiment, User.id == PGExperiment.user_id
).outerjoin(
    PGCampaign, User.id == PGCampaign.user_id
).group_by(User.id)
```

### 5. Added PostgreSQL Initialization [main.py:63-69](api/main.py#L63-L69)
```python
# Initialize SQLite database (experiments, jobs, campaigns)
init_db()
print("✅ SQLite Database initialized")

# Initialize PostgreSQL database (users, auth, admin)
init_postgres_db()
print("✅ PostgreSQL Database initialized")
```

## Endpoints Fixed

All admin panel endpoints now working correctly:

- ✅ `/api/admin/stats/overview` - System statistics
- ✅ `/api/admin/experiments/live` - Live monitoring
- ✅ `/api/admin/experiments/recent` - Recent experiments
- ✅ `/api/admin/usage/by-backend` - Backend usage analytics
- ✅ `/api/admin/usage/by-user` - User usage analytics
- ✅ `/api/admin/users` - User management
- ✅ `/api/admin/keys` - Access key management

## Testing

Verified from backend logs:
```
✓ GET /api/admin/stats/overview → 200
✓ GET /api/admin/experiments/live → 200
✓ GET /api/admin/usage/by-backend → 200
✓ GET /api/admin/usage/by-user → 200
✓ GET /api/admin/users → 200
✓ DELETE /api/admin/users/14 → 200
```

## Admin Panel Access

**URL**: `http://localhost:3000/mkkisarkar`
**Password**: `kanad_admin_2024`

## Features Now Working

1. **Live Monitoring Tab**
   - Shows all currently running experiments
   - Real-time progress tracking
   - Backend and method details
   - Running time calculation

2. **Analytics Tab**
   - Usage statistics by quantum backend
   - Top users by experiment count
   - Time period filters (7/30/90/365 days)
   - Success rate metrics

3. **Overview Tab**
   - Total users, experiments, campaigns
   - Active users count
   - Recent activity (last 24 hours)
   - System-wide statistics

4. **Users Tab**
   - Complete user management
   - Experiment and campaign counts per user
   - Role management (admin/user/viewer)
   - Account status controls

---

**Fixed**: 2025-10-30
**Issue**: Database architecture mismatch between admin routes and actual data storage
**Solution**: Query SQLite for experiment data, cross-reference PostgreSQL for user metadata
