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


def darken_color(rgb_tuple, factor=0.64):
    """Darken an RGB color by multiplying each channel by a factor."""
    r, g, b = rgb_tuple
    return (int(r * factor), int(g * factor), int(b * factor))
