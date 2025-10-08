"""BLADE artifact removal package."""

from .Artifact_remove import Artifact_remove

try:
    from .Artifact_remove_HD import Artifact_remove_HD, ArtifactRemoveHD
except ImportError:  # pragma: no cover
    Artifact_remove_HD = None
    ArtifactRemoveHD = None

__all__ = [
    name
    for name in ["Artifact_remove", "Artifact_remove_HD", "ArtifactRemoveHD"]
    if name in globals() and globals()[name] is not None
]
