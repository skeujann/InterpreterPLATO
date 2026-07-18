"""Pytest discovery wrappers for the standalone offline regression scripts."""

from __future__ import annotations

from InterpreterPLATO.scripts import test_compact_result_serialization
from InterpreterPLATO.scripts import test_expression_loading
from InterpreterPLATO.scripts import test_production_context


def test_expression_loading_regression() -> None:
    assert test_expression_loading.main() == 0


def test_production_context_regression() -> None:
    assert test_production_context.main() == 0


def test_compact_result_serialization_regression() -> None:
    assert test_compact_result_serialization.main() == 0
