"""Build a ChEMBL-target to HGNC-symbol mapping from the pinned ChEMBL r37 SQLite dump.

This is an identifier translation layer only.  It does not import ChEMBL
drug-target associations and must not be used as an additional target source.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


PACKAGE = Path(__file__).resolve().parents[1]
RAW_ROOT = PACKAGE / "data" / "raw"
PROCESSED_CONTEXT = PACKAGE / "data" / "processed" / "context"
HGNC_RELEASE = "2026-04-01"
CHEMBL_RELEASE = "r37"
DEFAULT_HGNC = RAW_ROOT / "gene_mapping" / "hgnc" / HGNC_RELEASE / f"hgnc_complete_set_{HGNC_RELEASE}.txt"
DEFAULT_CHEMBL_RAW = RAW_ROOT / "target_mapping" / "chembl" / CHEMBL_RELEASE
DEFAULT_SQLITE = DEFAULT_CHEMBL_RAW / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db"
DEFAULT_CHECKSUMS = DEFAULT_CHEMBL_RAW / "checksums.txt"
DEFAULT_ARCHIVE = DEFAULT_CHEMBL_RAW / "chembl_37_sqlite.tar.gz"
DEFAULT_OUTPUT = PROCESSED_CONTEXT / "chembl_r37_target_hgnc_mapping.tsv.gz"
DEFAULT_QC = PROCESSED_CONTEXT / "chembl_r37_target_hgnc_mapping_qc.json"
HASH_CHUNK_SIZE = 8 * 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _split_multi(value: Any) -> List[str]:
    text = str(value or "").strip()
    if not text:
        return []
    # HGNC complete-set fields commonly use pipe separators.  ChEMBL synonyms
    # are already one value per row, but accepting comma/semicolon is harmless.
    items: List[str] = []
    for part in text.replace("|", ";").replace(",", ";").split(";"):
        cleaned = " ".join(part.strip().split())
        if cleaned:
            items.append(cleaned)
    return list(dict.fromkeys(items))


def _norm_symbol(value: str) -> str:
    return str(value or "").strip().upper()


class HGNCResolver:
    """Resolve symbols and UniProt accessions into one fixed HGNC snapshot."""

    def __init__(self, path: Path):
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(f"HGNC file not found: {self.path}")
        self.sha256 = sha256_file(self.path)
        self.symbol_to_hgnc: Dict[str, str] = {}
        self.alias_to_symbols: Dict[str, set[str]] = defaultdict(set)
        self.uniprot_to_symbols: Dict[str, set[str]] = defaultdict(set)
        with self.path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"hgnc_id", "symbol", "status", "alias_symbol", "prev_symbol", "uniprot_ids"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"HGNC table missing columns: {sorted(missing)}")
            for row in reader:
                if str(row.get("status", "")).strip().lower() != "approved":
                    continue
                symbol = str(row.get("symbol", "")).strip()
                if not symbol:
                    continue
                normalized = _norm_symbol(symbol)
                self.symbol_to_hgnc[normalized] = str(row.get("hgnc_id", "")).strip()
                self.alias_to_symbols[normalized].add(normalized)
                for field in ("alias_symbol", "prev_symbol"):
                    for alias in _split_multi(row.get(field)):
                        self.alias_to_symbols[_norm_symbol(alias)].add(normalized)
                for accession in _split_multi(row.get("uniprot_ids")):
                    self.uniprot_to_symbols[accession.strip().upper()].add(normalized)
        if not self.symbol_to_hgnc:
            raise ValueError("HGNC table yielded no approved symbols")

    def resolve_symbol(self, value: str) -> Tuple[str, str, str]:
        normalized = _norm_symbol(value)
        if not normalized:
            return "", "", "empty_symbol"
        if normalized in self.symbol_to_hgnc:
            return normalized, self.symbol_to_hgnc[normalized], "approved_symbol"
        candidates = sorted(self.alias_to_symbols.get(normalized, set()))
        if len(candidates) == 1:
            symbol = candidates[0]
            return symbol, self.symbol_to_hgnc.get(symbol, ""), "unique_alias_or_previous_symbol"
        if len(candidates) > 1:
            return "", "", "ambiguous_hgnc_alias"
        return "", "", "symbol_not_in_hgnc"

    def resolve_uniprot(self, value: str) -> Tuple[str, str, str]:
        accession = str(value or "").strip().upper()
        if not accession:
            return "", "", "empty_uniprot"
        candidates = sorted(self.uniprot_to_symbols.get(accession, set()))
        if len(candidates) == 1:
            symbol = candidates[0]
            return symbol, self.symbol_to_hgnc.get(symbol, ""), "unique_uniprot"
        if len(candidates) > 1:
            return "", "", "ambiguous_uniprot"
        return "", "", "uniprot_not_in_hgnc"


def expected_archive_sha256(checksums_path: Path, filename: str) -> str:
    with Path(checksums_path).open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if str(row.get("file", "")).strip() == filename:
                return str(row.get("SHA256", "")).strip().lower()
    raise ValueError(f"{filename} not found in {checksums_path}")


def query_chembl_rows(sqlite_path: Path) -> Iterable[Dict[str, Any]]:
    query = """
        SELECT
            td.chembl_id AS target_chembl_id,
            td.pref_name AS target_name,
            td.target_type AS chembl_target_type,
            td.tax_id AS chembl_tax_id,
            td.organism AS chembl_organism,
            tc.component_id AS component_id,
            cs.accession AS uniprot_id,
            cs.component_type AS component_type,
            cs.tax_id AS component_tax_id,
            cs.organism AS component_organism,
            cs.description AS component_description,
            csyn.component_synonym AS component_synonym,
            csyn.syn_type AS syn_type
        FROM target_dictionary td
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
        LEFT JOIN component_synonyms csyn ON cs.component_id = csyn.component_id
        WHERE td.chembl_id IS NOT NULL
          AND (td.tax_id = 9606 OR cs.tax_id = 9606)
          AND (
                csyn.syn_type IN ('GENE_SYMBOL', 'GENE_SYMBOL_OTHER', 'UNIPROT')
                OR csyn.syn_type IS NULL
          )
        ORDER BY td.chembl_id, tc.component_id
    """
    if not sqlite_path.is_file():
        raise FileNotFoundError(f"ChEMBL SQLite DB not found: {sqlite_path}")
    connection = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        for row in connection.execute(query):
            yield dict(row)
    finally:
        connection.close()


def _choose_mapping(
    *,
    resolver: HGNCResolver,
    gene_symbols: Sequence[str],
    other_symbols: Sequence[str],
    uniprot_id: str,
) -> Dict[str, str]:
    attempts: List[Tuple[str, str, str, str]] = []
    for symbol in gene_symbols:
        canonical, hgnc_id, reason = resolver.resolve_symbol(symbol)
        attempts.append(("chembl_gene_symbol", symbol, canonical, reason))
        if canonical:
            return {
                "canonical_hgnc_symbol": canonical,
                "hgnc_id": hgnc_id,
                "source_gene_symbol": symbol,
                "source_gene_symbol_type": "GENE_SYMBOL",
                "mapping_method": "chembl_gene_symbol_to_hgnc",
                "mapping_status": "mapped",
                "mapping_reason": reason,
            }
    for symbol in other_symbols:
        canonical, hgnc_id, reason = resolver.resolve_symbol(symbol)
        attempts.append(("chembl_gene_symbol_other", symbol, canonical, reason))
        if canonical:
            return {
                "canonical_hgnc_symbol": canonical,
                "hgnc_id": hgnc_id,
                "source_gene_symbol": symbol,
                "source_gene_symbol_type": "GENE_SYMBOL_OTHER",
                "mapping_method": "chembl_gene_symbol_other_to_hgnc",
                "mapping_status": "mapped",
                "mapping_reason": reason,
            }
    canonical, hgnc_id, reason = resolver.resolve_uniprot(uniprot_id)
    attempts.append(("uniprot", uniprot_id, canonical, reason))
    if canonical:
        return {
            "canonical_hgnc_symbol": canonical,
            "hgnc_id": hgnc_id,
            "source_gene_symbol": "",
            "source_gene_symbol_type": "",
            "mapping_method": "uniprot_to_hgnc",
            "mapping_status": "mapped",
            "mapping_reason": reason,
        }
    reasons = [f"{kind}:{value}:{why}" for kind, value, _, why in attempts if value]
    return {
        "canonical_hgnc_symbol": "",
        "hgnc_id": "",
        "source_gene_symbol": gene_symbols[0] if gene_symbols else (other_symbols[0] if other_symbols else ""),
        "source_gene_symbol_type": "GENE_SYMBOL" if gene_symbols else ("GENE_SYMBOL_OTHER" if other_symbols else ""),
        "mapping_method": "unmapped",
        "mapping_status": "unmapped",
        "mapping_reason": "|".join(reasons) or "no_gene_symbol_or_uniprot_candidate",
    }


def build_mapping(sqlite_path: Path, hgnc_path: Path, output_path: Path, qc_path: Path) -> Dict[str, Any]:
    resolver = HGNCResolver(hgnc_path)
    grouped: Dict[Tuple[str, str], Dict[str, Any]] = {}
    target_component_counts: Dict[str, set[str]] = defaultdict(set)
    source_rows = 0
    for raw in query_chembl_rows(sqlite_path):
        source_rows += 1
        target_id = str(raw.get("target_chembl_id") or "").strip()
        component_id = str(raw.get("component_id") or "").strip()
        key = (target_id, component_id)
        record = grouped.setdefault(key, {
            "target_chembl_id": target_id,
            "target_name": str(raw.get("target_name") or "").strip(),
            "chembl_target_type": str(raw.get("chembl_target_type") or "").strip(),
            "chembl_organism": str(raw.get("chembl_organism") or "").strip(),
            "chembl_tax_id": str(raw.get("chembl_tax_id") or "").strip(),
            "component_id": component_id,
            "component_type": str(raw.get("component_type") or "").strip(),
            "component_uniprot_id": str(raw.get("uniprot_id") or "").strip(),
            "component_tax_id": str(raw.get("component_tax_id") or "").strip(),
            "component_organism": str(raw.get("component_organism") or "").strip(),
            "component_description": str(raw.get("component_description") or "").strip(),
            "gene_symbols": [],
            "gene_symbol_others": [],
        })
        if target_id and component_id:
            target_component_counts[target_id].add(component_id)
        synonym = str(raw.get("component_synonym") or "").strip()
        syn_type = str(raw.get("syn_type") or "").strip().upper()
        if not synonym:
            continue
        if syn_type == "GENE_SYMBOL" and synonym not in record["gene_symbols"]:
            record["gene_symbols"].append(synonym)
        elif syn_type == "GENE_SYMBOL_OTHER" and synonym not in record["gene_symbol_others"]:
            record["gene_symbol_others"].append(synonym)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "target_chembl_id", "target_name", "chembl_target_type", "chembl_organism", "chembl_tax_id",
        "component_id", "component_type", "component_uniprot_id", "component_tax_id", "component_organism",
        "component_description", "source_gene_symbol", "source_gene_symbol_type",
        "canonical_hgnc_symbol", "hgnc_id", "mapping_method", "mapping_status", "mapping_reason",
        "target_component_count",
    ]
    counters: Counter[str] = Counter()
    unique_targets: set[str] = set()
    mapped_targets: set[str] = set()
    mapped_single_component_targets: set[str] = set()
    temporary_output = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    try:
        with temporary_output.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0,
            ) as compressed_handle:
                with io.TextIOWrapper(
                    compressed_handle, encoding="utf-8", newline="",
                ) as handle:
                    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
                    writer.writeheader()
                    for _, record in sorted(grouped.items()):
                        mapping = _choose_mapping(
                            resolver=resolver,
                            gene_symbols=record.pop("gene_symbols"),
                            other_symbols=record.pop("gene_symbol_others"),
                            uniprot_id=str(record.get("component_uniprot_id", "")),
                        )
                        target_id = str(record["target_chembl_id"])
                        component_count = len(target_component_counts.get(target_id, set()))
                        row = {
                            **record,
                            **mapping,
                            "target_component_count": component_count,
                        }
                        writer.writerow(row)
                        unique_targets.add(target_id)
                        counters[f"status:{mapping['mapping_status']}"] += 1
                        counters[f"method:{mapping['mapping_method']}"] += 1
                        counters[f"target_type:{record['chembl_target_type']}"] += 1
                        if mapping["mapping_status"] == "mapped":
                            mapped_targets.add(target_id)
                            if component_count == 1:
                                mapped_single_component_targets.add(target_id)
        os.replace(str(temporary_output), str(output_path))
    finally:
        if temporary_output.exists():
            temporary_output.unlink()

    output_sha256 = sha256_file(output_path)
    qc = {
        "schema_version": "1.0",
        "source": "ChEMBL",
        "release": CHEMBL_RELEASE,
        "chembl_sqlite_path": str(sqlite_path.resolve()),
        "chembl_sqlite_size_bytes": sqlite_path.stat().st_size,
        "chembl_sqlite_sha256": sha256_file(sqlite_path),
        "hgnc_release": HGNC_RELEASE,
        "hgnc_path": str(Path(hgnc_path).resolve()),
        "hgnc_sha256": resolver.sha256,
        "output_path": str(output_path.resolve()),
        "output_size_bytes": output_path.stat().st_size,
        "output_sha256": output_sha256,
        "built_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_sql_rows": source_rows,
        "component_rows": len(grouped),
        "unique_targets": len(unique_targets),
        "mapped_targets": len(mapped_targets),
        "mapped_single_component_targets": len(mapped_single_component_targets),
        "counters": dict(counters),
        "semantics": (
            "Identifier translation only: ChEMBL target IDs from PLATO are mapped "
            "to HGNC symbols; no ChEMBL drug-target associations are imported."
        ),
    }
    qc_path.write_text(json.dumps(qc, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return qc


def main() -> int:
    parser = argparse.ArgumentParser(description="Build ChEMBL r37 target ID to HGNC mapping")
    parser.add_argument("--sqlite", type=Path, default=DEFAULT_SQLITE)
    parser.add_argument("--hgnc", type=Path, default=DEFAULT_HGNC)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--checksums", type=Path, default=DEFAULT_CHECKSUMS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qc", type=Path, default=DEFAULT_QC)
    parser.add_argument("--skip-archive-check", action="store_true")
    args = parser.parse_args()

    if not args.skip_archive_check:
        expected = expected_archive_sha256(args.checksums, args.archive.name)
        observed = sha256_file(args.archive)
        if observed.lower() != expected.lower():
            raise SystemExit(
                f"ChEMBL archive SHA-256 mismatch: expected {expected}, observed {observed}"
            )

    qc = build_mapping(args.sqlite, args.hgnc, args.output, args.qc)
    print(json.dumps(qc, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
