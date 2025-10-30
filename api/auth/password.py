"""
Password Hashing and Verification for Kanad Platform
Uses bcrypt for secure password hashing
"""

import bcrypt
from typing import Optional
import re


def hash_password(password: str) -> str:
    """
    Hash a password using bcrypt

    Args:
        password: Plain text password

    Returns:
        Hashed password string
    """
    # Generate salt and hash password
    salt = bcrypt.gensalt(rounds=12)
    hashed = bcrypt.hashpw(password.encode("utf-8"), salt)
    return hashed.decode("utf-8")


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """
    Verify a password against its hash

    Args:
        plain_password: Plain text password to verify
        hashed_password: Hashed password from database

    Returns:
        True if password matches, False otherwise
    """
    try:
        return bcrypt.checkpw(
            plain_password.encode("utf-8"), hashed_password.encode("utf-8")
        )
    except Exception:
        return False


def validate_password_strength(password: str) -> tuple[bool, Optional[str]]:
    """
    Validate password meets security requirements

    Requirements:
    - Minimum 8 characters
    - At least one uppercase letter
    - At least one lowercase letter
    - At least one digit
    - At least one special character

    Args:
        password: Password to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if len(password) < 8:
        return False, "Password must be at least 8 characters long"

    if not re.search(r"[A-Z]", password):
        return False, "Password must contain at least one uppercase letter"

    if not re.search(r"[a-z]", password):
        return False, "Password must contain at least one lowercase letter"

    if not re.search(r"\d", password):
        return False, "Password must contain at least one digit"

    if not re.search(r"[!@#$%^&*(),.?\":{}|<>]", password):
        return False, "Password must contain at least one special character (!@#$%^&*(),.?\":{}|<>)"

    return True, None


def generate_temporary_password(length: int = 12) -> str:
    """
    Generate a random temporary password

    Args:
        length: Length of password (default 12)

    Returns:
        Random password string meeting security requirements
    """
    import secrets
    import string

    # Ensure password has required character types
    uppercase = secrets.choice(string.ascii_uppercase)
    lowercase = secrets.choice(string.ascii_lowercase)
    digit = secrets.choice(string.digits)
    special = secrets.choice("!@#$%^&*()")

    # Fill remaining length with random characters
    all_chars = string.ascii_letters + string.digits + "!@#$%^&*()"
    remaining = "".join(secrets.choice(all_chars) for _ in range(length - 4))

    # Combine and shuffle
    password_list = list(uppercase + lowercase + digit + special + remaining)
    secrets.SystemRandom().shuffle(password_list)

    return "".join(password_list)


def is_password_compromised(password: str) -> bool:
    """
    Check if password is in common password list (basic check)

    This is a simple implementation. For production, consider using
    HaveIBeenPwned API or a more comprehensive password list.

    Args:
        password: Password to check

    Returns:
        True if password is common/compromised, False otherwise
    """
    common_passwords = {
        "password",
        "password123",
        "123456",
        "12345678",
        "qwerty",
        "abc123",
        "monkey",
        "1234567",
        "letmein",
        "trustno1",
        "dragon",
        "baseball",
        "111111",
        "iloveyou",
        "master",
        "sunshine",
        "ashley",
        "bailey",
        "passw0rd",
        "shadow",
        "123123",
        "654321",
        "superman",
        "qazwsx",
        "michael",
        "football",
    }

    return password.lower() in common_passwords


def validate_and_hash_password(password: str) -> tuple[Optional[str], Optional[str]]:
    """
    Validate password strength and hash if valid

    Args:
        password: Plain text password

    Returns:
        Tuple of (hashed_password, error_message)
        If valid: (hashed_password, None)
        If invalid: (None, error_message)
    """
    # Check password strength
    is_valid, error = validate_password_strength(password)
    if not is_valid:
        return None, error

    # Check if password is compromised
    if is_password_compromised(password):
        return None, "This password is too common. Please choose a more secure password."

    # Hash and return
    hashed = hash_password(password)
    return hashed, None


class PasswordValidator:
    """
    Password validation class with configurable rules
    """

    def __init__(
        self,
        min_length: int = 8,
        max_length: int = 128,
        require_uppercase: bool = True,
        require_lowercase: bool = True,
        require_digit: bool = True,
        require_special: bool = True,
        check_compromised: bool = True,
    ):
        self.min_length = min_length
        self.max_length = max_length
        self.require_uppercase = require_uppercase
        self.require_lowercase = require_lowercase
        self.require_digit = require_digit
        self.require_special = require_special
        self.check_compromised = check_compromised

    def validate(self, password: str) -> tuple[bool, Optional[str]]:
        """
        Validate password against configured rules

        Args:
            password: Password to validate

        Returns:
            Tuple of (is_valid, error_message)
        """
        # Check length
        if len(password) < self.min_length:
            return False, f"Password must be at least {self.min_length} characters long"

        if len(password) > self.max_length:
            return False, f"Password must be no more than {self.max_length} characters long"

        # Check character requirements
        if self.require_uppercase and not re.search(r"[A-Z]", password):
            return False, "Password must contain at least one uppercase letter"

        if self.require_lowercase and not re.search(r"[a-z]", password):
            return False, "Password must contain at least one lowercase letter"

        if self.require_digit and not re.search(r"\d", password):
            return False, "Password must contain at least one digit"

        if self.require_special and not re.search(r"[!@#$%^&*(),.?\":{}|<>]", password):
            return (
                False,
                "Password must contain at least one special character (!@#$%^&*(),.?\":{}|<>)",
            )

        # Check compromised
        if self.check_compromised and is_password_compromised(password):
            return False, "This password is too common. Please choose a more secure password."

        return True, None


# Default validator instance
default_validator = PasswordValidator()
