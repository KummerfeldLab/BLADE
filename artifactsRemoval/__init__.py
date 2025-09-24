"""
Package initializer for artifactsRemoval.
Avoid hard imports that are optional or absent in some environments to keep
submodule imports (like Artifact_remove_HD) working without requiring others.
"""

# Export non-HD Artifact_remove if available, but don't fail if not.
try:
    from .Artifact_remove import Artifact_remove  # noqa: F401
except Exception:
    pass

# Export load_data if available, but don't fail if not (e.g., missing helper).
try:
    from .load_sample import load_data  # noqa: F401
except Exception:
    pass
