"""Pytest-native wrappers for the original strict contract assertions."""

from __future__ import annotations

from pathlib import Path

import pytest

import InterpreterPLATO.targets as target_module
from InterpreterPLATO.scripts import test_compact as compact


PLATO_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data" / "raw" / "plato" / "r37" / "golden_fixture"
    / "gemcitabine_target_fishing_response_2026-07-16.json"
)


def test_no_production_schema_fixture_is_approved() -> None:
    assert target_module.VERIFIED_PLATO_SCHEMA_FINGERPRINTS == ()


def test_target_parsing_filters() -> None:
    compact.test_target_parsing_filters()


@pytest.mark.skipif(
    not PLATO_FIXTURE.is_file(),
    reason="external PLATO fixture is intentionally excluded from the code-only repository",
)
def test_real_plato_r37_fixture_contract() -> None:
    compact.test_real_plato_r37_fixture_contract()


@pytest.mark.skipif(
    not PLATO_FIXTURE.is_file(),
    reason="external PLATO fixture is intentionally excluded from the code-only repository",
)
def test_exact_candidate_is_rejected_even_when_raw_hash_is_trusted(tmp_path: Path) -> None:
    compact.test_exact_candidate_is_rejected_even_when_raw_hash_is_trusted(tmp_path)


@pytest.mark.skipif(
    not PLATO_FIXTURE.is_file(),
    reason="external PLATO fixture is intentionally excluded from the code-only repository",
)
def test_stale_schema_grant_is_revoked_when_production_allowlist_is_empty(tmp_path: Path) -> None:
    compact.test_stale_schema_grant_is_revoked_when_production_allowlist_is_empty(tmp_path)


def test_strict_api_rejects_local_and_mutable_modes(tmp_path: Path) -> None:
    compact.test_strict_api_rejects_local_and_mutable_modes(tmp_path)


def test_unverified_and_untrusted_plato_exports_fail(tmp_path: Path) -> None:
    compact.test_unverified_and_untrusted_plato_exports_fail(tmp_path)


def test_strict_positive_fixture_and_integrity(tmp_path: Path) -> None:
    compact.test_strict_positive_fixture_and_integrity(tmp_path)


def test_invalid_smiles_fail_closed(tmp_path: Path) -> None:
    compact.test_invalid_smiles_fail_closed(tmp_path)
