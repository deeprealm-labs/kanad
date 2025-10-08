"""
Secure credential management with encryption.

Uses Fernet (symmetric encryption) to encrypt/decrypt cloud provider credentials.
"""

from cryptography.fernet import Fernet
from api.config import settings
import base64
import logging

logger = logging.getLogger(__name__)


class CredentialsManager:
    """
    Manage encrypted storage of cloud provider credentials.

    Uses Fernet symmetric encryption (AES-128 in CBC mode).
    """

    def __init__(self):
        """Initialize with encryption key from settings."""
        # Encryption key must be 32 url-safe base64-encoded bytes
        key = settings.ENCRYPTION_KEY.encode()

        # Validate key format
        try:
            self.cipher = Fernet(key)
        except Exception as e:
            logger.error(f"Invalid encryption key: {e}")
            raise ValueError(
                "ENCRYPTION_KEY must be a valid Fernet key. "
                "Generate with: python -c 'from cryptography.fernet import Fernet; print(Fernet.generate_key().decode())'"
            )

    def encrypt_token(self, token: str) -> str:
        """
        Encrypt an API token.

        Args:
            token: Plain text token

        Returns:
            Base64-encoded encrypted token
        """
        if not token:
            return ""

        try:
            encrypted = self.cipher.encrypt(token.encode())
            return base64.b64encode(encrypted).decode()
        except Exception as e:
            logger.error(f"Encryption failed: {e}")
            raise

    def decrypt_token(self, encrypted_token: str) -> str:
        """
        Decrypt an API token.

        Args:
            encrypted_token: Base64-encoded encrypted token

        Returns:
            Plain text token
        """
        if not encrypted_token:
            return ""

        try:
            encrypted = base64.b64decode(encrypted_token.encode())
            decrypted = self.cipher.decrypt(encrypted)
            return decrypted.decode()
        except Exception as e:
            logger.error(f"Decryption failed: {e}")
            raise

    def encrypt_credentials(self, credentials: dict) -> dict:
        """
        Encrypt all credential fields.

        Args:
            credentials: Dict with plaintext credentials

        Returns:
            Dict with encrypted credentials
        """
        encrypted = {}
        for key, value in credentials.items():
            if value:
                encrypted[f"encrypted_{key}"] = self.encrypt_token(value)
            else:
                encrypted[f"encrypted_{key}"] = ""

        return encrypted

    def decrypt_credentials(self, encrypted_credentials: dict) -> dict:
        """
        Decrypt all credential fields.

        Args:
            encrypted_credentials: Dict with encrypted credentials

        Returns:
            Dict with plaintext credentials
        """
        decrypted = {}
        for key, value in encrypted_credentials.items():
            if key.startswith("encrypted_") and value:
                original_key = key.replace("encrypted_", "")
                decrypted[original_key] = self.decrypt_token(value)
            elif not key.startswith("encrypted_"):
                # Pass through non-encrypted fields
                decrypted[key] = value

        return decrypted


# Global instance
credentials_manager = CredentialsManager()
