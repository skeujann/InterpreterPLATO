"""Regression tests for vectorized long/wide expression-file loading."""

from __future__ import annotations

import gzip
from hashlib import sha256
from pathlib import Path
import sys
import tempfile


sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO.context import CancerContext


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _assert_clean_audit(context: CancerContext) -> None:
    audit = context.audit_runtime_resources()
    assert audit["mutated_resources"] == []
    assert audit["findings"] == []


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="interpreterplato-expression-") as temporary:
        root = Path(temporary)

        # Wide gene x cell-line format, including duplicate genes and invalid
        # numeric values. Valid duplicates retain the last observation.
        wide = root / "wide.csv.gz"
        with gzip.open(wide, "wt", encoding="utf-8", newline="") as handle:
            handle.write(
                "gene,A549,MCF7\n"
                "egfr,1,2\n"
                " SRC ,3,not-a-number\n"
                "EGFR,4,inf\n"
                "BRAF,-inf,5\n"
            )
        wide_context = CancerContext().load_expression(
            wide,
            source_name="synthetic-wide",
            version="test-v1",
            background_complete=True,
            detection_threshold=1.0,
            taxid="9606",
            identifier_namespace="hgnc_symbol",
            mapping_version="HGNC-test",
        )
        assert wide_context.expression == {
            "A549": {"EGFR": 4.0, "SRC": 3.0},
            "MCF7": {"EGFR": 2.0, "BRAF": 5.0},
        }
        assert all(
            type(value) is float
            for profile in wide_context.expression.values()
            for value in profile.values()
        )
        wide_metadata = wide_context.resource_metadata["expression"]
        assert wide_metadata["cell_lines"] == ["A549", "MCF7"]
        assert wide_metadata["gene_counts"] == {"A549": 2, "MCF7": 2}
        assert set(wide_metadata["artifact_sha256_by_cell_line"].values()) == {_sha256(wide)}
        assert wide_metadata["background_complete_by_cell_line"] == {
            "A549": True, "MCF7": True,
        }
        assert wide_metadata["detection_threshold_by_cell_line"] == {
            "A549": 1.0, "MCF7": 1.0,
        }
        _assert_clean_audit(wide_context)

        # Long format retains the same invalid-value filtering, last-row-wins
        # behavior, and the historical distinction between an incoming empty
        # line and a line with at least one valid observation.
        long = root / "long.csv"
        long.write_text(
            "gene,expression,cell_line\n"
            "egfr,1,A549\n"
            "EGFR,4,A549\n"
            "SRC,not-a-number,A549\n"
            "BRAF,3,MCF7\n"
            "MAPK1,inf,MCF7\n"
            "AKT1,bad,EMPTY\n",
            encoding="utf-8",
        )
        long_context = CancerContext().load_expression(
            long,
            source_name="synthetic-long",
            version="test-v2",
            background_complete=False,
            taxid="9606",
            identifier_namespace="hgnc_symbol",
            mapping_version="HGNC-test",
        )
        assert long_context.expression == {
            "A549": {"EGFR": 4.0},
            "MCF7": {"BRAF": 3.0},
            "EMPTY": {},
        }
        long_metadata = long_context.resource_metadata["expression"]
        assert long_metadata["cell_lines"] == ["A549", "EMPTY", "MCF7"]
        assert long_metadata["gene_counts"] == {"A549": 1, "MCF7": 1, "EMPTY": 0}
        assert set(long_metadata["source_by_cell_line"]) == {"A549", "MCF7"}
        assert set(long_metadata["artifact_sha256_by_cell_line"].values()) == {_sha256(long)}
        _assert_clean_audit(long_context)

        no_line = root / "long_without_cell_line.csv"
        no_line.write_text("gene,value\nakt1,1.5\nAKT1,2.5\n", encoding="utf-8")
        no_line_context = CancerContext().load_expression(no_line, cell_line="CUSTOM")
        assert no_line_context.expression == {"CUSTOM": {"AKT1": 2.5}}
        _assert_clean_audit(no_line_context)

    print("expression loading tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
