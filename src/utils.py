"""Shared utility functions for SQANTI-browser."""

import re


def normalize_trix_token(value):
    """Normalize strings for Trix search tokens (lowercase, underscores)."""
    if value is None:
        return ""
    token = str(value).strip().lower()
    if not token:
        return ""
    if token == '+':
        return 'plus'
    if token == '-':
        return 'minus'
    token = re.sub(r'[^\w]+', '_', token)
    token = re.sub(r'_+', '_', token)
    token = token.strip('_')
    return token
