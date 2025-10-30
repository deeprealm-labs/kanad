"""
Email OTP (One-Time Password) Service for Kanad Platform
Handles email verification using 6-digit OTP codes
"""

import os
import secrets
import string
from datetime import datetime, timedelta
from typing import Optional, Dict, Any
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from sqlalchemy.orm import Session

# Email configuration
SMTP_HOST = os.getenv("SMTP_HOST", "smtp.gmail.com")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER", "")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
SMTP_FROM_EMAIL = os.getenv("SMTP_FROM_EMAIL", SMTP_USER)
SMTP_FROM_NAME = os.getenv("SMTP_FROM_NAME", "Kanad Platform")

# OTP configuration
OTP_LENGTH = 6
OTP_EXPIRY_MINUTES = 10
OTP_MAX_ATTEMPTS = 3


def generate_otp(length: int = OTP_LENGTH) -> str:
    """
    Generate a random numeric OTP

    Args:
        length: Length of OTP (default 6)

    Returns:
        Random numeric OTP string
    """
    return "".join(secrets.choice(string.digits) for _ in range(length))


def create_otp(db: Session, email: str) -> tuple[str, datetime]:
    """
    Create and store OTP in database

    Args:
        db: Database session
        email: Email address to send OTP to

    Returns:
        Tuple of (otp, expires_at)
    """
    from api.core.database_postgres import EmailVerification

    # Generate OTP
    otp = generate_otp()
    expires_at = datetime.utcnow() + timedelta(minutes=OTP_EXPIRY_MINUTES)

    # Store in database
    verification = EmailVerification(
        email=email, otp=otp, expires_at=expires_at, is_used=False
    )

    db.add(verification)
    db.commit()

    return otp, expires_at


def verify_otp(db: Session, email: str, otp: str) -> tuple[bool, Optional[str]]:
    """
    Verify OTP for email

    Args:
        db: Database session
        email: Email address
        otp: OTP code to verify

    Returns:
        Tuple of (is_valid, error_message)
    """
    from api.core.database_postgres import EmailVerification

    # Find most recent unused OTP for this email
    verification = (
        db.query(EmailVerification)
        .filter(
            EmailVerification.email == email,
            EmailVerification.otp == otp,
            EmailVerification.is_used == False,
        )
        .order_by(EmailVerification.created_at.desc())
        .first()
    )

    if not verification:
        return False, "Invalid OTP code"

    # Check if expired
    if datetime.utcnow() > verification.expires_at:
        return False, "OTP code has expired"

    # Mark as used
    verification.is_used = True
    db.commit()

    return True, None


def send_otp_email(email: str, otp: str, user_name: Optional[str] = None) -> bool:
    """
    Send OTP via email

    Args:
        email: Recipient email address
        otp: OTP code to send
        user_name: Optional user name for personalization

    Returns:
        True if email sent successfully, False otherwise
    """
    try:
        # Create message
        msg = MIMEMultipart("alternative")
        msg["Subject"] = "Verify Your Email - Kanad Platform"
        msg["From"] = f"{SMTP_FROM_NAME} <{SMTP_FROM_EMAIL}>"
        msg["To"] = email

        # Create HTML and plain text versions
        greeting = f"Hello {user_name}," if user_name else "Hello,"

        text_body = f"""
{greeting}

Thank you for registering with Kanad Platform!

Your verification code is: {otp}

This code will expire in {OTP_EXPIRY_MINUTES} minutes.

If you didn't request this code, please ignore this email.

Best regards,
Kanad Platform Team
        """

        html_body = f"""
<!DOCTYPE html>
<html>
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0;">
  </head>
  <body style="margin: 0; padding: 0; font-family: Georgia, 'Times New Roman', serif; background: #f5f5f5;">
    <table width="100%" cellpadding="0" cellspacing="0" border="0" style="background: #f5f5f5;">
      <tr>
        <td align="center" style="padding: 20px 0;">
          <table width="600" cellpadding="0" cellspacing="0" border="0" style="background: #ffffff; max-width: 600px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
            <!-- Header -->
            <tr>
              <td style="background: #000000; padding: 50px 30px; text-align: center;">
                <h1 style="color: #D2691E; font-size: 64px; margin: 0; font-weight: 700; letter-spacing: 3px; font-family: Georgia, 'Times New Roman', serif;">kanad</h1>
                <p style="color: #ffffff; font-size: 13px; margin: 15px 0 0 0; letter-spacing: 3px; text-transform: lowercase; opacity: 0.7; font-family: Georgia, serif;">quantum chemistry platform</p>
              </td>
            </tr>

            <!-- Content -->
            <tr>
              <td style="padding: 40px 35px; background: #ffffff;">
                <p style="font-size: 20px; color: #1a1a1a; margin: 0 0 20px 0; font-family: Georgia, serif;">{greeting}</p>
                <p style="font-size: 15px; color: #6b7280; margin: 0 0 15px 0; line-height: 1.7; font-family: Georgia, serif;">
                  Thank you for joining Kanad. To complete your registration and unlock quantum chemistry simulations,
                  please verify your email address using the code below.
                </p>

                <!-- OTP -->
                <table width="100%" cellpadding="0" cellspacing="0" border="0" style="margin: 35px 0;">
                  <tr>
                    <td align="center">
                      <p style="font-size: 11px; color: #6b7280; text-transform: uppercase; letter-spacing: 1.5px; margin: 0 0 15px 0; font-family: Arial, sans-serif; font-weight: 600;">YOUR VERIFICATION CODE</p>
                      <table cellpadding="0" cellspacing="0" border="0" style="background: #000000; border: 2px solid #D2691E;">
                        <tr>
                          <td style="padding: 30px 40px;">
                            <p style="font-size: 48px; font-weight: 700; color: #D2691E; letter-spacing: 14px; font-family: 'Courier New', Courier, monospace; margin: 0;">{otp}</p>
                          </td>
                        </tr>
                      </table>
                    </td>
                  </tr>
                </table>

                <!-- Warning -->
                <table width="100%" cellpadding="0" cellspacing="0" border="0" style="margin: 25px 0;">
                  <tr>
                    <td style="background: #fff8f0; border: 2px solid #D2691E; padding: 18px; font-size: 14px; color: #1a1a1a; font-family: Georgia, serif;">
                      <p style="margin: 0;"><strong style="color: #D2691E;">⏱ Time Sensitive:</strong> This code will expire in <strong>{OTP_EXPIRY_MINUTES} minutes</strong>.</p>
                    </td>
                  </tr>
                </table>

                <!-- Security -->
                <table width="100%" cellpadding="0" cellspacing="0" border="0" style="margin: 25px 0;">
                  <tr>
                    <td style="background: #f9fafb; border-left: 4px solid #000000; padding: 18px; font-size: 13px; color: #6b7280; line-height: 1.6; font-family: Georgia, serif;">
                      <p style="margin: 0;"><strong style="color: #1a1a1a;">🔒 Security Notice:</strong> If you didn't request this code, please ignore this email. Never share your verification code with anyone.</p>
                    </td>
                  </tr>
                </table>

                <p style="font-size: 15px; color: #6b7280; margin: 30px 0 0 0; line-height: 1.7; font-family: Georgia, serif;">
                  Best regards,<br><strong style="color: #1a1a1a;">The Kanad Team</strong>
                </p>
              </td>
            </tr>

            <!-- Footer -->
            <tr>
              <td style="background: #000000; color: #6b7280; padding: 30px; text-align: center;">
                <p style="color: #D2691E; font-size: 32px; margin: 0 0 15px 0; font-weight: 700; letter-spacing: 2px; font-family: Georgia, serif;">kanad</p>
                <p style="font-size: 12px; color: #6b7280; margin: 0; line-height: 1.6; font-family: Georgia, serif;">
                  © 2024 Kanad. All rights reserved.<br>
                  This is an automated message. Please do not reply to this email.
                </p>
              </td>
            </tr>
          </table>
        </td>
      </tr>
    </table>
  </body>
</html>
        """

        # Attach both versions
        part1 = MIMEText(text_body, "plain")
        part2 = MIMEText(html_body, "html")
        msg.attach(part1)
        msg.attach(part2)

        # Send email
        if not SMTP_USER or not SMTP_PASSWORD:
            print("⚠️  SMTP credentials not configured. OTP email not sent.")
            print(f"OTP for {email}: {otp}")
            return True  # Return True for development

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            server.starttls()
            server.login(SMTP_USER, SMTP_PASSWORD)
            server.send_message(msg)

        print(f"✓ OTP email sent to {email}")
        return True

    except Exception as e:
        print(f"✗ Failed to send OTP email: {str(e)}")
        return False


def resend_otp(db: Session, email: str) -> tuple[bool, Optional[str]]:
    """
    Resend OTP to email

    Args:
        db: Database session
        email: Email address

    Returns:
        Tuple of (success, error_message)
    """
    from api.core.database_postgres import EmailVerification

    # Check recent OTP attempts
    recent_attempts = (
        db.query(EmailVerification)
        .filter(
            EmailVerification.email == email,
            EmailVerification.created_at
            > datetime.utcnow() - timedelta(minutes=5),
        )
        .count()
    )

    if recent_attempts >= 3:
        return False, "Too many OTP requests. Please wait 5 minutes before trying again."

    # Create and send new OTP
    otp, expires_at = create_otp(db, email)
    success = send_otp_email(email, otp)

    if not success:
        return False, "Failed to send OTP email. Please try again later."

    return True, None


def cleanup_expired_otps(db: Session) -> int:
    """
    Clean up expired OTP records from database

    Args:
        db: Database session

    Returns:
        Number of records deleted
    """
    from api.core.database_postgres import EmailVerification

    # Delete OTPs older than 24 hours
    cutoff_time = datetime.utcnow() - timedelta(hours=24)

    deleted = (
        db.query(EmailVerification)
        .filter(EmailVerification.created_at < cutoff_time)
        .delete()
    )

    db.commit()

    print(f"✓ Cleaned up {deleted} expired OTP records")
    return deleted


def get_otp_status(db: Session, email: str) -> Dict[str, Any]:
    """
    Get OTP verification status for email

    Args:
        db: Database session
        email: Email address

    Returns:
        Dictionary with OTP status information
    """
    from api.core.database_postgres import EmailVerification

    # Get most recent OTP
    latest_otp = (
        db.query(EmailVerification)
        .filter(EmailVerification.email == email)
        .order_by(EmailVerification.created_at.desc())
        .first()
    )

    if not latest_otp:
        return {
            "has_otp": False,
            "is_expired": False,
            "is_used": False,
            "expires_at": None,
            "created_at": None,
        }

    is_expired = datetime.utcnow() > latest_otp.expires_at

    return {
        "has_otp": True,
        "is_expired": is_expired,
        "is_used": latest_otp.is_used,
        "expires_at": latest_otp.expires_at.isoformat(),
        "created_at": latest_otp.created_at.isoformat(),
    }


# Email templates for different scenarios
def send_welcome_email(email: str, user_name: str) -> bool:
    """
    Send welcome email after successful registration

    Args:
        email: User's email address
        user_name: User's name

    Returns:
        True if sent successfully
    """
    try:
        msg = MIMEMultipart("alternative")
        msg["Subject"] = "Welcome to Kanad Platform!"
        msg["From"] = f"{SMTP_FROM_NAME} <{SMTP_FROM_EMAIL}>"
        msg["To"] = email

        html_body = f"""
<html>
  <body style="font-family: Arial, sans-serif;">
    <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
      <h2 style="color: #ff6b35;">Welcome to Kanad Platform, {user_name}!</h2>
      <p>Your account has been successfully created and verified.</p>
      <p>You can now start exploring quantum chemistry simulations and analysis.</p>
      <p>If you have any questions, feel free to reach out to our support team.</p>
      <p>Best regards,<br>Kanad Platform Team</p>
    </div>
  </body>
</html>
        """

        msg.attach(MIMEText(html_body, "html"))

        if SMTP_USER and SMTP_PASSWORD:
            with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
                server.starttls()
                server.login(SMTP_USER, SMTP_PASSWORD)
                server.send_message(msg)

        return True

    except Exception as e:
        print(f"✗ Failed to send welcome email: {str(e)}")
        return False
