"""BLADE artifact removal package."""

from .Artifact_remove import Artifact_remove

try:
    from .Artifact_remove_HD import Artifact_remove_HD, ArtifactRemoveHD
except ImportError:  # pragma: no cover
    Artifact_remove_HD = None
    ArtifactRemoveHD = None

try:
    from .iterative_edge_cleanup import run_iterative_edge_cleanup
except ImportError:  # pragma: no cover
    run_iterative_edge_cleanup = None

__all__ = [
    name
    for name in ["Artifact_remove", "Artifact_remove_HD", "ArtifactRemoveHD", "run_iterative_edge_cleanup"]
    if name in globals() and globals()[name] is not None
]
