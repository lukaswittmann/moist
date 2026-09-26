"""Python API for moist solvation models."""

import re

from .density import GaussianBasis, InternalDensity
from .library import get_version_string, get_banner

from .interface import (
    AnchorGradient,
    Cavity,
    CavityDROP,
    CavityField,
    CavityISwiG,
    CavitySnapshot,
    CavitySnapshotDROP,
    Coupling,
    CouplingRequest,
    DensityResponse,
    PotentialAdjointResponse,
    GaussianMomentRequest,
    GostshypAmplitudeResponse,
    PointPotentialRequest,
    GaussianPotentialRequest,
    Response,
    IsodensitySource,
    ModelComponentCOSMO,
    ModelComponentCPCM,
    ModelComponentGOSTSHYP,
    ModelComponentPV,
    PCMSolver,
    SolvationModel,
    SolvationModelComponent,
    Structure,
)

from .parameters import (
    CFCParameters,
    DROPParameters,
    ISwiGParameters,
    IsodensityParameters,
    ModelParameters,
    PCMParameters,
    SvdWParameters,
)
from .radii import BondiRadii, COSMORadii, CPCMRadii, CustomRadii, D3Radii, Radii, SMDRadii
from .configuration import CFC, DROP, ISwiG, Isodensity, SvdW, CavityConfiguration, LevelSet

__all__ = [
    "GaussianBasis", "InternalDensity", "get_version_string", "get_banner",
    "DROPParameters",
    "ISwiGParameters",
    "SvdWParameters",
    "CFCParameters",
    "IsodensityParameters",
    "PCMParameters",
    "ModelParameters",
    "Radii",
    "CPCMRadii",
    "SMDRadii",
    "D3Radii",
    "COSMORadii",
    "BondiRadii",
    "CustomRadii",
    "LevelSet",
    "CavityConfiguration",
    "DROP",
    "ISwiG",
    "SvdW",
    "CFC",
    "Isodensity",
    "AnchorGradient",
    "Cavity",
    "CavityDROP",
    "CavityField",
    "CavityISwiG",
    "CavitySnapshot",
    "CavitySnapshotDROP",
    "Coupling",
    "CouplingRequest",
    "DensityResponse",
    "PotentialAdjointResponse",
    "GaussianMomentRequest",
    "GostshypAmplitudeResponse",
    "PointPotentialRequest",
    "GaussianPotentialRequest",
    "Response",
    "IsodensitySource",
    "ModelComponentCOSMO",
    "ModelComponentCPCM",
    "ModelComponentGOSTSHYP",
    "ModelComponentPV",
    "PCMSolver",
    "SolvationModel",
    "SolvationModelComponent",
    "Structure",
]


def _pep440(version: str) -> str:
    """Map the native ``1.0.0-alpha.1`` style tag to PEP 440 (``1.0.0a1``)."""
    match = re.fullmatch(r"(\d+\.\d+\.\d+)(?:-(alpha|beta|rc)\.?(\d+))?", version)
    if match is None:
        return version
    base, stage, number = match.groups()
    if stage is None:
        return base
    return base + {"alpha": "a", "beta": "b", "rc": "rc"}[stage] + number


#: Version of the linked native library, the single source of truth.
__version__ = _pep440(get_version_string())
