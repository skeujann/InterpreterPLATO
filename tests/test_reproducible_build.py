"""Byte-level reproducibility checks for generated runtime artifacts."""

from __future__ import annotations

from pathlib import Path

from InterpreterPLATO.scripts import build_chembl_target_mapping as builder


def test_chembl_mapping_gzip_is_reproducible(
    tmp_path: Path, monkeypatch,
) -> None:
    hgnc = tmp_path / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tstatus\talias_symbol\tprev_symbol\tuniprot_ids\n"
        "HGNC:3236\tEGFR\tApproved\tERBB\t\tP00533\n",
        encoding="utf-8",
    )
    sqlite = tmp_path / "chembl.db"
    sqlite.write_bytes(b"synthetic sqlite identity")
    rows = [{
        "target_chembl_id": "CHEMBL203",
        "target_name": "Epidermal growth factor receptor",
        "chembl_target_type": "SINGLE PROTEIN",
        "chembl_tax_id": 9606,
        "chembl_organism": "Homo sapiens",
        "component_id": 1,
        "uniprot_id": "P00533",
        "component_type": "PROTEIN",
        "component_tax_id": 9606,
        "component_organism": "Homo sapiens",
        "component_description": "EGFR",
        "component_synonym": "EGFR",
        "syn_type": "GENE_SYMBOL",
    }]
    monkeypatch.setattr(builder, "query_chembl_rows", lambda _: list(rows))

    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    builder.build_mapping(sqlite, hgnc, first, tmp_path / "first.qc.json")
    builder.build_mapping(sqlite, hgnc, second, tmp_path / "second.qc.json")

    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes()[4:8] == b"\x00\x00\x00\x00"
