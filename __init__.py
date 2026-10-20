"""Compact PLATO-aware synergy mechanism interpreter."""

__version__ = "5.5.2"

from .advanced import AdvancedResult, explain_advanced
from .context import CancerContext, ContextScore
from .interpreter import (
    CombinationExplanation,
    CounterfactualHypothesisEngine,
    DrugTargetExplanation,
    ExplanationPipeline,
    ModuleExplanation,
    ModuleHit,
    SynergyInterpreter,
)
from .mechanism import MechanismResult, classify_mechanism
from .targets import (
    PlatoPredictionCache,
    PlatoTargetEvidence,
    inspect_plato_payload,
    parse_plato_json,
    standardize_smiles,
)

__all__ = [
    "SynergyInterpreter",
    "ExplanationPipeline",
    "CombinationExplanation",
    "DrugTargetExplanation",
    "ModuleHit",
    "ModuleExplanation",
    "CounterfactualHypothesisEngine",
    "PlatoPredictionCache",
    "PlatoTargetEvidence",
    "parse_plato_json",
    "inspect_plato_payload",
    "standardize_smiles",
    "CancerContext",
    "ContextScore",
    "MechanismResult",
    "classify_mechanism",
    "AdvancedResult",
    "explain_advanced",
]
