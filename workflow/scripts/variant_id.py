import hashlib


MAX_VARIANT_ID_LENGTH = 200


def safe_variant_id(variant_id):
    """Return a stable ID that is safe to use as a profile filename."""
    if len(variant_id) <= MAX_VARIANT_ID_LENGTH:
        return variant_id
    return "sv_" + hashlib.sha256(variant_id.encode()).hexdigest()[:20]
