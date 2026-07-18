"""Shared pytest isolation for the historical regression suite."""

from __future__ import annotations

from pathlib import Path
import sys

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PACKAGE_ROOT.parent
SCRIPTS_ROOT = PACKAGE_ROOT / "scripts"
for candidate in (WORKSPACE_ROOT, SCRIPTS_ROOT):
    value = str(candidate)
    if value not in sys.path:
        sys.path.insert(0, value)

import InterpreterPLATO.targets as target_module  # noqa: E402


@pytest.fixture(autouse=True)
def restore_plato_schema_allowlist() -> None:
    """Do not let one fixture's process-local schema grant leak to another test."""
    original = tuple(target_module.VERIFIED_PLATO_SCHEMA_FINGERPRINTS)
    try:
        yield
    finally:
        target_module.VERIFIED_PLATO_SCHEMA_FINGERPRINTS = original
