# Google OAuth Setup Guide for Kanad

## Error: "Access blocked: Authorization Error" - origin_mismatch

This error occurs when the redirect URI in your application doesn't match what's configured in Google Cloud Console.

---

## Quick Fix

### Step 1: Check Your Current Settings

Your backend currently uses:
- **Redirect URI**: `http://localhost:3000/auth/google/callback`
- **Client ID**: Set via `GOOGLE_CLIENT_ID` environment variable
- **Client Secret**: Set via `GOOGLE_CLIENT_SECRET` environment variable

### Step 2: Configure Google Cloud Console

1. **Go to Google Cloud Console**
   - Visit: https://console.cloud.google.com/

2. **Select or Create a Project**
   - Click project dropdown (top left)
   - Create new project: "Kanad Platform"

3. **Enable Google+ API**
   - Go to: APIs & Services → Library
   - Search for "Google+ API"
   - Click "Enable"

4. **Configure OAuth Consent Screen**
   - Go to: APIs & Services → OAuth consent screen
   - Choose "External" (unless you have Google Workspace)
   - Fill in required fields:
     - **App name**: Kanad Quantum Chemistry Platform
     - **User support email**: your-email@gmail.com
     - **Developer contact**: your-email@gmail.com
   - Click "Save and Continue"
   - Skip scopes (default is fine)
   - Add test users: your-email@gmail.com
   - Click "Save and Continue"

5. **Create OAuth Credentials**
   - Go to: APIs & Services → Credentials
   - Click "Create Credentials" → "OAuth client ID"
   - Application type: "Web application"
   - Name: "Kanad Web Client"

   **Authorized JavaScript origins**:
   Add ALL origins you'll use:
   ```
   http://localhost:3000
   http://127.0.0.1:3000
   http://172.171.222.16
   https://your-app.vercel.app
   ```

   **Authorized redirect URIs**:
   Add ALL redirect URIs:
   ```
   http://localhost:3000/auth/google/callback
   http://127.0.0.1:3000/auth/google/callback
   http://172.171.222.16/auth/google/callback
   https://your-app.vercel.app/auth/google/callback
   ```

   - Click "Create"
   - **IMPORTANT**: Copy the Client ID and Client Secret

---

## Step 3: Update Backend Environment

### Local Development

Edit `.env` file:
```bash
cd /home/mk/deeprealm/kanad
nano .env
```

Add these lines (replace with your actual values):
```env
# Google OAuth Configuration
GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your-client-secret
GOOGLE_REDIRECT_URI=http://localhost:3000/auth/google/callback
```

### Azure Production

SSH into Azure VM:
```bash
ssh kanadmin@172.171.222.16
cd /opt/kanad
nano .env
```

Add these lines:
```env
# Google OAuth Configuration
GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your-client-secret
GOOGLE_REDIRECT_URI=http://172.171.222.16/auth/google/callback
```

Restart backend:
```bash
sudo systemctl restart kanad
```

---

## Step 4: Update Frontend Configuration

The frontend needs to know the Google Client ID for client-side authentication.

### Local Development

Edit `web/.env.local`:
```bash
cd /home/mk/deeprealm/kanad/web
nano .env.local
```

Add:
```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api
NEXT_PUBLIC_GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
```

### Production (Vercel)

In Vercel Dashboard → Settings → Environment Variables:
```
NEXT_PUBLIC_API_URL=http://172.171.222.16/api
NEXT_PUBLIC_GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
```

---

## Step 5: Frontend Google OAuth Component

Check if your frontend has Google OAuth button. If not, here's how it should look:

### Example: Login page with Google OAuth

```tsx
import { useEffect } from 'react';

declare global {
  interface Window {
    google: any;
  }
}

export default function LoginPage() {
  useEffect(() => {
    // Load Google Sign-In script
    const script = document.createElement('script');
    script.src = 'https://accounts.google.com/gsi/client';
    script.async = true;
    script.defer = true;
    document.body.appendChild(script);

    script.onload = () => {
      window.google.accounts.id.initialize({
        client_id: process.env.NEXT_PUBLIC_GOOGLE_CLIENT_ID,
        callback: handleGoogleCallback,
      });

      window.google.accounts.id.renderButton(
        document.getElementById('google-signin-button'),
        { theme: 'outline', size: 'large', text: 'signin_with' }
      );
    };

    return () => {
      document.body.removeChild(script);
    };
  }, []);

  const handleGoogleCallback = async (response: any) => {
    try {
      // Send ID token to backend
      const result = await fetch('http://localhost:8000/api/auth/google', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ id_token: response.credential }),
      });

      const data = await result.json();

      if (result.ok) {
        // Store tokens
        localStorage.setItem('kanad_access_token', data.access_token);
        localStorage.setItem('kanad_refresh_token', data.refresh_token);

        // Redirect to dashboard
        window.location.href = '/dashboard';
      } else {
        alert(data.detail || 'Google login failed');
      }
    } catch (error) {
      console.error('Google login error:', error);
      alert('Failed to login with Google');
    }
  };

  return (
    <div>
      <h1>Login to Kanad</h1>
      <div id="google-signin-button"></div>
    </div>
  );
}
```

---

## Testing Google OAuth

### 1. Check Backend Configuration
```bash
curl http://localhost:8000/api/auth/status
```

Should return OAuth configuration status.

### 2. Test Google Login Flow

1. Start backend:
   ```bash
   cd /home/mk/deeprealm/kanad
   source env/bin/activate
   python -m uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
   ```

2. Start frontend:
   ```bash
   cd web
   npm run dev
   ```

3. Open browser: http://localhost:3000

4. Click "Sign in with Google"

5. Should redirect to Google OAuth consent screen

6. After approving, should redirect back with token

### 3. Verify Authentication

Check if user was created:
```bash
# Connect to database
psql -U postgres -d kanad_db

# Query users
SELECT id, email, google_id, is_verified, created_at FROM users;

# Exit
\q
```

---

## Common Issues & Solutions

### Issue 1: "origin_mismatch" Error

**Symptom**:
```
Error 400: origin_mismatch
The redirect URI in the request doesn't match the ones authorized for the OAuth client.
```

**Solution**:
1. Check the URL in your browser address bar
2. Go to Google Cloud Console → Credentials
3. Edit your OAuth client
4. Add the EXACT origin to "Authorized JavaScript origins"
5. Add the EXACT redirect URI to "Authorized redirect URIs"

**Example**:
If you're at `http://localhost:3000`, add:
- Origin: `http://localhost:3000`
- Redirect: `http://localhost:3000/auth/google/callback`

### Issue 2: "Access blocked: This app's request is invalid"

**Symptom**: Google shows "This app hasn't been verified" warning

**Solution**:
1. Go to OAuth consent screen in Google Cloud Console
2. Add your email as a test user
3. For production, submit app for verification (can take 1-2 weeks)
4. Alternatively, use "Continue" button during testing

### Issue 3: "redirect_uri_mismatch"

**Symptom**: Error mentions redirect_uri mismatch

**Solution**:
Check that `GOOGLE_REDIRECT_URI` in backend `.env` matches exactly what's in Google Cloud Console.

Common mistakes:
- Trailing slash: `http://localhost:3000/` vs `http://localhost:3000`
- HTTP vs HTTPS
- Port number missing
- Path incorrect

### Issue 4: "invalid_client"

**Symptom**: Backend logs show "invalid_client" error

**Solution**:
1. Verify `GOOGLE_CLIENT_ID` and `GOOGLE_CLIENT_SECRET` are correct
2. Check for extra spaces or newlines in `.env` file
3. Ensure credentials match the OAuth client in Google Cloud Console

### Issue 5: "Access Key Required"

**Symptom**: "Access key required for new Google user"

**Solution**:
For new users signing up with Google, you need an access key:

```bash
# Create access key (as admin)
curl -X POST http://localhost:8000/api/admin/keys \
  -H "Authorization: Bearer YOUR_ADMIN_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "description": "Google OAuth access key",
    "max_uses": 100,
    "expires_in_days": 365
  }'
```

Then provide the access key during Google OAuth:
```javascript
const result = await fetch('http://localhost:8000/api/auth/google', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    id_token: response.credential,
    access_key: 'your-access-key-here'
  }),
});
```

---

## Environment Variables Reference

### Backend (.env)

```env
# Required for Google OAuth
GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your-client-secret
GOOGLE_REDIRECT_URI=http://localhost:3000/auth/google/callback

# For production (Azure)
# GOOGLE_REDIRECT_URI=http://172.171.222.16/auth/google/callback
# For Vercel
# GOOGLE_REDIRECT_URI=https://your-app.vercel.app/auth/google/callback
```

### Frontend (.env.local)

```env
# Required for client-side Google Sign-In
NEXT_PUBLIC_GOOGLE_CLIENT_ID=your-client-id.apps.googleusercontent.com
NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

---

## Security Best Practices

1. **Never commit credentials**
   - Add `.env` to `.gitignore`
   - Use environment variables on production

2. **Use HTTPS in production**
   - Set up SSL/TLS on Azure
   - Update redirect URIs to use HTTPS

3. **Limit origins**
   - Only add origins you control
   - Remove test origins in production

4. **Rotate secrets regularly**
   - Generate new client secret every 90 days
   - Update in all environments

5. **Monitor OAuth usage**
   - Check Google Cloud Console → APIs & Services → Metrics
   - Watch for suspicious activity

---

## Production Deployment Checklist

- [ ] Google Cloud Console project created
- [ ] OAuth consent screen configured
- [ ] Production origins added:
  - [ ] Vercel frontend URL
  - [ ] Azure backend IP/domain
- [ ] Redirect URIs added:
  - [ ] Vercel callback URL
  - [ ] Azure callback URL
- [ ] Backend environment variables set (Azure)
- [ ] Frontend environment variables set (Vercel)
- [ ] SSL/HTTPS enabled on backend
- [ ] Test Google login on production
- [ ] Verify user creation in database
- [ ] Monitor for errors

---

## Multi-Environment Setup

### Development
```env
GOOGLE_REDIRECT_URI=http://localhost:3000/auth/google/callback
```

### Staging (if applicable)
```env
GOOGLE_REDIRECT_URI=https://staging.yourdomain.com/auth/google/callback
```

### Production
```env
GOOGLE_REDIRECT_URI=https://kanad.yourdomain.com/auth/google/callback
```

**Note**: You can use the SAME Client ID and Secret across all environments by adding all redirect URIs to the same OAuth client in Google Cloud Console.

---

## Debugging Steps

1. **Check backend logs**:
   ```bash
   # Local
   # Check terminal where uvicorn is running

   # Azure
   ssh kanadmin@172.171.222.16
   sudo journalctl -u kanad -f
   ```

2. **Check frontend console**:
   - Open browser DevTools (F12)
   - Go to Console tab
   - Look for Google OAuth errors

3. **Test backend OAuth endpoint**:
   ```bash
   curl -X POST http://localhost:8000/api/auth/google \
     -H "Content-Type: application/json" \
     -d '{"id_token":"test-token"}'
   ```

4. **Verify environment variables**:
   ```bash
   # Backend
   cd /home/mk/deeprealm/kanad
   grep GOOGLE .env

   # Frontend
   cd web
   grep GOOGLE .env.local
   ```

---

## Additional Resources

- **Google OAuth Documentation**: https://developers.google.com/identity/protocols/oauth2
- **Google Sign-In for Websites**: https://developers.google.com/identity/sign-in/web
- **OAuth 2.0 Playground**: https://developers.google.com/oauthplayground/

---

## Summary

✅ **To fix the "origin_mismatch" error**:

1. Go to https://console.cloud.google.com/apis/credentials
2. Edit your OAuth 2.0 Client ID
3. Add your current browser URL origin to "Authorized JavaScript origins"
4. Add the callback URL to "Authorized redirect URIs"
5. Save and try again

**Common fix**:
If you're testing locally at `http://localhost:3000`:
- Add origin: `http://localhost:3000`
- Add redirect: `http://localhost:3000/auth/google/callback`

**For Azure production**:
- Add origin: `http://172.171.222.16`
- Add redirect: `http://172.171.222.16/auth/google/callback`

**For Vercel production**:
- Add origin: `https://your-app.vercel.app`
- Add redirect: `https://your-app.vercel.app/auth/google/callback`

---

*Last Updated: 2025-10-30*
*Status: Troubleshooting Guide*
