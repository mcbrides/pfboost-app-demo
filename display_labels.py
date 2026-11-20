"""Short browsing labels; full construct boundaries remain in the source metadata."""
import re


def display_variant(label):
    value = str(label or '').strip()
    # Remove only plain sequence intervals/lengths, never mutation positions or
    # circular-permutation names (PM23, P54–55), background residues, or WT*.
    value = re.sub(r'(?<![A-Za-z0-9])\d+\s*[–-]\s*\d+(?![A-Za-z0-9])', '', value)
    value = re.sub(r'\b\d+\s*-?\s*aa\b', '', value, flags=re.I)
    value = re.sub(r'\bfull-length\s+(?=[A-Z]\d+[A-Z])', '', value)
    value = re.sub(r'\bfull-length\b', '', value)
    # Clean punctuation left by deleting an interval inside a qualification.
    value = re.sub(r'\(\s*[;,]?\s*', '(', value)
    value = re.sub(r'\s*[;,]?\s*\)', ')', value)
    value = re.sub(r'\(\s*\)', '', value)
    value = re.sub(r'\s+', ' ', value).strip()
    if value == 'Mature':
        return 'Mature core'
    if value in {'', 'core'}:
        return 'Core construct' if label else ''
    value = re.sub(r'\(core\)', '', value).strip()
    value = re.sub(r'\(mature\)', '(mature core)', value)
    value = re.sub(r'\bMature\s+(?=\+|$)', 'Mature core ', value).strip()
    return value
