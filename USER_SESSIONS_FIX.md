# User Sessions API Fix

## Issue
The `/api/users/me/sessions` endpoint was returning 500 Internal Server Error due to Pydantic validation failure.

## Root Cause
**Response Model Mismatch**: The `SessionInfo` Pydantic model expected a field called `last_activity`, but the endpoint was returning `last_used` from the database session object.

### Error Message
```
fastapi.exceptions.ResponseValidationError: 1 validation errors:
  {'type': 'missing', 'loc': ('response', 0, 'last_activity'), 'msg': 'Field required', 'input': {'id': 21, 'created_at': datetime.datetime(...), 'expires_at': datetime.datetime(...), 'last_used': datetime.datetime(...), 'is_current': True}}
```

The error showed that the response had `last_used` but the model required `last_activity`.

## Fix Applied

### 1. Updated SessionInfo Model in [users.py:56-61](api/routes/users.py#L56-L61)
```python
class SessionInfo(BaseModel):
    id: int
    created_at: datetime
    expires_at: datetime
    last_activity: Optional[datetime]  # Changed from last_used to match frontend expectations
    is_current: bool
```

### 2. Updated get_user_sessions Response in [users.py:225-234](api/routes/users.py#L225-L234)
```python
return [
    {
        "id": session.id,
        "created_at": session.created_at,
        "expires_at": session.expires_at,
        "last_activity": session.last_used,  # Map last_used to last_activity for frontend
        "is_current": session.id == current_session_id,
    }
    for session in sessions
]
```

## Verification

After the fix was applied and the server reloaded automatically, the endpoint started working correctly:

```
✓ GET /api/users/me/sessions → 200
INFO:     127.0.0.1:35554 - "GET /api/users/me/sessions HTTP/1.1" 200 OK
```

## Related Endpoints

All user profile endpoints are now working correctly:
- ✅ `/api/users/me` - Get current user profile
- ✅ `/api/users/me/sessions` - Get user sessions
- ✅ `/api/users/me/statistics` - Get user statistics
- ✅ `/api/users/change-password` - Change password
- ✅ `/api/users/me` (DELETE) - Delete account

---

**Fixed**: 2025-10-30
**Issue**: Pydantic model field name mismatch between database field and API response
**Solution**: Map `last_used` database field to `last_activity` API response field
