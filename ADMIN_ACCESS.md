# Admin Panel Access

## Security Fixes Applied

### Critical Security Vulnerabilities Fixed:

1. **Dashboard Routes Not Protected** ✅ FIXED
   - Previously: Anyone could access `/dashboard` routes without authentication
   - Now: All `/dashboard` routes require login authentication
   - Implementation: Added auth check in `dashboard/layout.tsx`

2. **Admin Panel Exposed in Sidebar** ✅ FIXED
   - Previously: Admin link visible in sidebar (even though backend was protected)
   - Now: Admin panel completely separated to obscure URL
   - Implementation: Removed from sidebar, moved to independent route

3. **Obscure Admin URL** ✅ IMPLEMENTED
   - Admin panel not accessible at predictable `/admin` URL
   - Moved to non-obvious path for security through obscurity
   - This prevents automated scanning and casual discovery

## Admin Panel Access

### URL:
```
http://localhost:3000/mkkisarkar
```

⚠️ **KEEP THIS URL SECRET** - Do not share publicly or commit to version control

### Password:
```
kanad_admin_2024
```

### Features:
- System statistics overview
- Access key management (create, deactivate, view usage)
- User management (view, edit roles, deactivate accounts)
- Live experiment monitoring (real-time updates)
- Usage analytics (by backend, by user, timelines)

## Security Measures Implemented:

1. **Frontend Password Protection**
   - Password stored in session storage (valid for browser session)
   - Password required on every new browser session
   - Clean UI with Kanad branding

2. **Backend API Protection**
   - All admin API endpoints require `require_admin` dependency
   - JWT token validation
   - Role-based access control (admin role required)

3. **Dashboard Authentication**
   - Login required for all `/dashboard` routes
   - Automatic redirect to home page if not authenticated
   - Loading states during auth checks

## Changing the Admin Password

To change the admin password, edit the file:
```
web/src/app/admin/page.tsx
```

Line 105:
```typescript
const ADMIN_PASSWORD = "kanad_admin_2024";
```

Change `"kanad_admin_2024"` to your desired password.

## Production Recommendations:

For production deployment, consider:

1. **Environment Variable Password**
   ```typescript
   const ADMIN_PASSWORD = process.env.NEXT_PUBLIC_ADMIN_PASSWORD;
   ```

2. **Hash-Based Authentication**
   - Store password hash instead of plain text
   - Use bcrypt or similar for verification

3. **Backend Admin Login Endpoint**
   - Create dedicated admin login API endpoint
   - Store admin sessions in database
   - Implement rate limiting for login attempts

4. **Two-Factor Authentication**
   - Add 2FA for admin access
   - Use TOTP (Google Authenticator, Authy)

5. **IP Whitelist**
   - Restrict admin panel access to specific IPs
   - Use firewall rules or middleware

6. **Audit Logging**
   - Log all admin actions
   - Track who accessed what and when
   - Alert on suspicious activities

## Testing the Fixes:

1. **Test Dashboard Protection:**
   - Open browser in incognito mode
   - Navigate to `http://localhost:3000/dashboard`
   - Should redirect to home page
   - Login with your account
   - Should now show dashboard

2. **Test Admin Panel:**
   - Navigate to `http://localhost:3000/admin`
   - Should show password prompt
   - Enter password: `kanad_admin_2024`
   - Should show admin dashboard with all features

3. **Test Session Persistence:**
   - Access admin panel and enter password
   - Navigate away and back to `/admin`
   - Should stay authenticated (session storage)
   - Close browser and reopen
   - Should require password again

## Admin Panel Features:

### Overview Tab
- Total users, experiments, campaigns statistics
- Active users count
- Success rates
- Quick action buttons

### Access Keys Tab
- Create new registration keys
- Set max uses and expiration
- View usage progress
- Deactivate or delete keys

### Users Tab
- View all registered users
- Filter by role, verified, active status
- Edit user roles and status
- View user activity stats
- Delete user accounts

### Live Monitoring Tab
- Real-time experiment tracking
- Auto-refresh every 5 seconds
- Progress bars and runtime
- Backend and method details

### Analytics Tab
- Usage by quantum backend
- Top users by activity
- Success rate metrics
- Time period filters (7/30/90/365 days)

---

**Last Updated:** 2025-10-30
**Security Level:** Medium (suitable for development, needs hardening for production)
