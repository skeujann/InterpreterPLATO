"""Build versioned production context artifacts from local raw data.

The builder is intentionally independent from the runtime interpreter.  It
normalizes every gene identifier through the fixed HGNC 2026-04-01 snapshot,
records every unresolved/ambiguous decision, and writes formats already
accepted by :class:`InterpreterPLATO.context.CancerContext`.

Examples
--------
    python InterpreterPLATO/scripts/build_context_artifacts.py hgnc
    python InterpreterPLATO/scripts/build_context_artifacts.py ppi
    python InterpreterPLATO/scripts/build_context_artifacts.py depmap
    python InterpreterPLATO/scripts/build_context_artifacts.py reactome
    python InterpreterPLATO/scripts/build_context_artifacts.py cosmic
    python InterpreterPLATO/scripts/build_context_artifacts.py all

No network access is used.  Production HGNC loading is checksum-locked.  The
``verify_hgnc_checksum=False`` function argument exists only so unit tests can
use a small synthetic HGNC fixture; the CLI never disables verification.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, TextIO, Tuple

import networkx as nx


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = PACKAGE_ROOT / "data" / "raw"
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "data" / "processed" / "context"

HGNC_RELEASE = "2026-04-01"
HGNC_EXPECTED_SHA256 = "ad19966c08ccebd7b0fad25d553a5db8770a188fef9cf93fc24dd74796c45d34"
STRING_RELEASE = "12.0"
DEPMAP_RELEASE = "26Q1"
REACTOME_RELEASE = "96"
COSMIC_RELEASE = "v101"

DEFAULT_HGNC = RAW_ROOT / "gene_mapping" / "hgnc" / HGNC_RELEASE / f"hgnc_complete_set_{HGNC_RELEASE}.txt"
DEFAULT_STRING_INFO = RAW_ROOT / "ppi" / "string" / "9606.protein.info.v12.0.txt"
DEFAULT_STRING_LINKS = RAW_ROOT / "ppi" / "string" / "9606.protein.links.full.v12.0.txt"
DEFAULT_DEPMAP_EXPRESSION = (
    RAW_ROOT / "cell_line_context" / "depmap" / DEPMAP_RELEASE
    / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)
DEFAULT_DEPMAP_MODEL = RAW_ROOT / "cell_line_context" / "depmap" / DEPMAP_RELEASE / "Model.csv"
DEFAULT_DEPMAP_PROFILES = RAW_ROOT / "cell_line_context" / "depmap" / DEPMAP_RELEASE / "OmicsProfiles.csv"
DEFAULT_SCORE_COMBOS = PACKAGE_ROOT / "data" / "scoreri_cell_lines.tsv"
DEFAULT_REACTOME = RAW_ROOT / "pathways" / "reactome" / "ReactomePathways.gmt"
DEFAULT_COSMIC = (
    RAW_ROOT / "disease" / "cosmic" / COSMIC_RELEASE
    / "Cosmic_CancerGeneCensus_v101_GRCh38.tsv"
)

PPI_OUTPUT_NAME = f"string_v{STRING_RELEASE}_hgnc_{HGNC_RELEASE}.tsv.gz"
EXPRESSION_OUTPUT_NAME = f"depmap_{DEPMAP_RELEASE}_scoreri_log2tpm.csv.gz"
REACTOME_OUTPUT_NAME = f"reactome_v{REACTOME_RELEASE}_hgnc_{HGNC_RELEASE}.gmt"
COSMIC_OUTPUT_NAME = f"cosmic_cgc_{COSMIC_RELEASE}_hgnc_{HGNC_RELEASE}.tsv"
COSMIC_TIER1_OUTPUT_NAME = f"cosmic_cgc_{COSMIC_RELEASE}_pan_cancer_tier1.txt"


# These are deliberate, reviewable decisions for the four ScoreRI names that
# do not uniquely match DepMap after punctuation/case normalization.  SR is
# explicitly blocked: the NCI-60 SR line is not the SR-786 line in DepMap.
EXPLICIT_SCORE_CELL_LINE_ALIASES: Dict[str, Dict[str, str]] = {
    "786-0": {
        "model_id": "ACH-000649",
        "expected_depmap_name": "786-O",
        "reason": "ScoreRI uses digit zero; explicitly mapped to DepMap 786-O",
    },
    "OVCAR3": {
        "model_id": "ACH-000001",
        "expected_depmap_name": "NIH:OVCAR-3",
        "reason": "ScoreRI shorthand explicitly mapped to NIH:OVCAR-3",
    },
    "U251": {
        "model_id": "ACH-000232",
        "expected_depmap_name": "U-251 MG",
        "reason": "explicitly mapped to parental U-251 MG, not U-251 MG DM",
    },
}
EXPLICIT_SCORE_CELL_LINE_REJECTIONS: Dict[str, str] = {
    "SR": (
        "unresolved by design: NCI-60 SR must not be guessed as DepMap SR-786; "
        "an authoritative RRID/crosswalk is required"
    ),
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_file(path: Path, label: str) -> Path:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def _file_record(path: Path) -> Dict[str, Any]:
    path = _require_file(path, "input")
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "sha256": _sha256(path),
    }


def _assert_unchanged(before: Mapping[str, Any]) -> None:
    path = Path(str(before["path"]))
    after_size = path.stat().st_size
    after_hash = _sha256(path)
    if after_size != before["size_bytes"] or after_hash != before["sha256"]:
        raise RuntimeError(f"input changed while being processed: {path}")


def _tmp_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


@contextmanager
def _atomic_text(path: Path, newline: str = "") -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline=newline) as handle:
            yield handle
        os.replace(str(temporary), str(path))
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def _atomic_gzip_text(path: Path) -> Iterator[TextIO]:
    """Write deterministic gzip (mtime=0) and atomically replace destination."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _tmp_path(path)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                    yield handle
        os.replace(str(temporary), str(path))
    finally:
        if temporary.exists():
            temporary.unlink()


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    with _atomic_text(path, newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    opener = _atomic_gzip_text if path.suffix.lower() == ".gz" else _atomic_text
    with opener(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _split_hgnc_multi(value: Any) -> List[str]:
    return [part.strip().upper() for part in str(value or "").split("|") if part.strip()]


def _clean_entrez(value: Any) -> str:
    text = str(value or "").strip()
    if re.fullmatch(r"\d+\.0", text):
        text = text[:-2]
    return text if text.isdigit() else ""


@dataclass(frozen=True)
class GeneMappingDecision:
    source_symbol: str
    source_entrez_id: str
    canonical_symbol: str
    hgnc_id: str
    status: str
    method: str
    candidates: Tuple[str, ...] = ()
    reason: str = ""

    def row(self) -> Dict[str, Any]:
        result = asdict(self)
        result["candidates"] = "|".join(self.candidates)
        return result


class HGNCMapper:
    """Strict mapper over one checksum-locked HGNC complete-set snapshot."""

    def __init__(self, path: Path, expected_sha256: Optional[str] = HGNC_EXPECTED_SHA256):
        self.path = _require_file(Path(path), "HGNC mapping")
        self.sha256 = _sha256(self.path)
        if expected_sha256 and self.sha256.lower() != expected_sha256.lower():
            raise ValueError(
                f"HGNC {HGNC_RELEASE} checksum mismatch: expected {expected_sha256}, got {self.sha256}"
            )
        if HGNC_RELEASE not in self.path.name and expected_sha256:
            raise ValueError(f"HGNC filename does not declare fixed release {HGNC_RELEASE}: {self.path.name}")

        self.approved: Dict[str, str] = {}
        self.hgnc_ids: Dict[str, str] = {}
        self.entrez: Dict[str, Set[str]] = defaultdict(set)
        self.aliases: Dict[str, Set[str]] = defaultdict(set)
        self.nonapproved_rows = 0
        required = {"hgnc_id", "symbol", "status", "entrez_id", "alias_symbol", "prev_symbol"}
        with self.path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"HGNC table missing columns: {sorted(missing)}")
            for row_number, row in enumerate(reader, start=2):
                if str(row.get("status", "")).strip().lower() != "approved":
                    self.nonapproved_rows += 1
                    continue
                symbol = str(row.get("symbol", "")).strip().upper()
                if not symbol:
                    raise ValueError(f"HGNC approved row {row_number} has no symbol")
                if symbol in self.approved:
                    raise ValueError(f"duplicate approved HGNC symbol at row {row_number}: {symbol}")
                self.approved[symbol] = symbol
                self.hgnc_ids[symbol] = str(row.get("hgnc_id", "")).strip()
                entrez = _clean_entrez(row.get("entrez_id"))
                if entrez:
                    self.entrez[entrez].add(symbol)
                for field in ("alias_symbol", "prev_symbol"):
                    for alias in _split_hgnc_multi(row.get(field)):
                        self.aliases[alias].add(symbol)
        if not self.approved:
            raise ValueError("HGNC table yielded no approved symbols")

    def map(self, symbol: Any = "", entrez_id: Any = "") -> GeneMappingDecision:
        raw_symbol = str(symbol or "").strip()
        normalized_symbol = raw_symbol.upper()
        normalized_entrez = _clean_entrez(entrez_id)

        if normalized_entrez:
            entrez_candidates = tuple(sorted(self.entrez.get(normalized_entrez, set())))
            if len(entrez_candidates) > 1:
                return GeneMappingDecision(
                    raw_symbol, normalized_entrez, "", "", "rejected", "ambiguous_entrez",
                    entrez_candidates, "Entrez ID maps to multiple approved HGNC symbols",
                )
            if len(entrez_candidates) == 1:
                canonical = entrez_candidates[0]
                reason = ""
                method = "entrez_unique"
                if normalized_symbol and normalized_symbol in self.approved and normalized_symbol != canonical:
                    method = "entrez_unique_over_symbol_conflict"
                    reason = f"Entrez-preferred mapping conflicts with approved source symbol {normalized_symbol}"
                elif normalized_symbol:
                    alias_candidates = self.aliases.get(normalized_symbol, set())
                    if len(alias_candidates) == 1 and canonical not in alias_candidates:
                        method = "entrez_unique_over_alias_conflict"
                        reason = "Entrez-preferred mapping conflicts with unique alias/previous-symbol mapping"
                return GeneMappingDecision(
                    raw_symbol, normalized_entrez, canonical, self.hgnc_ids.get(canonical, ""),
                    "mapped", method, (canonical,), reason,
                )

        if normalized_symbol in self.approved:
            reason = "Entrez ID not present in HGNC snapshot; exact approved symbol used" if normalized_entrez else ""
            return GeneMappingDecision(
                raw_symbol, normalized_entrez, normalized_symbol, self.hgnc_ids.get(normalized_symbol, ""),
                "mapped", "approved_exact", (normalized_symbol,), reason,
            )
        alias_candidates = tuple(sorted(self.aliases.get(normalized_symbol, set())))
        if len(alias_candidates) == 1:
            canonical = alias_candidates[0]
            reason = "Entrez ID not present in HGNC snapshot; unique alias/previous symbol used" if normalized_entrez else ""
            return GeneMappingDecision(
                raw_symbol, normalized_entrez, canonical, self.hgnc_ids.get(canonical, ""),
                "mapped", "unique_alias_or_previous", alias_candidates, reason,
            )
        if len(alias_candidates) > 1:
            return GeneMappingDecision(
                raw_symbol, normalized_entrez, "", "", "rejected", "ambiguous_alias_or_previous",
                alias_candidates, "identifier maps to multiple approved HGNC symbols",
            )
        return GeneMappingDecision(
            raw_symbol, normalized_entrez, "", "", "rejected", "unmapped", (),
            "identifier is not an approved symbol or unique HGNC alias/previous symbol",
        )

    def qc(self) -> Dict[str, Any]:
        return {
            "release": HGNC_RELEASE,
            "path": str(self.path.resolve()),
            "sha256": self.sha256,
            "approved_symbol_count": len(self.approved),
            "unique_entrez_id_count": sum(len(values) == 1 for values in self.entrez.values()),
            "ambiguous_entrez_id_count": sum(len(values) > 1 for values in self.entrez.values()),
            "alias_or_previous_key_count": len(self.aliases),
            "ambiguous_alias_or_previous_key_count": sum(len(values) > 1 for values in self.aliases.values()),
            "nonapproved_rows_excluded": self.nonapproved_rows,
            "mapping_policy": [
                "unique Entrez mapping first",
                "approved symbol exact match second",
                "unique alias/previous-symbol match third",
                "ambiguous identifiers rejected without fallback guessing",
            ],
        }


def build_hgnc_reports(mapper: HGNCMapper, output_dir: Path) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    ambiguity_path = output_dir / "hgnc_ambiguous_identifiers.tsv.gz"
    rows: List[Dict[str, Any]] = []
    for value, candidates in sorted(mapper.entrez.items()):
        if len(candidates) > 1:
            rows.append({"identifier_type": "entrez_id", "identifier": value, "candidates": "|".join(sorted(candidates))})
    for value, candidates in sorted(mapper.aliases.items()):
        if len(candidates) > 1:
            rows.append({"identifier_type": "alias_or_previous", "identifier": value, "candidates": "|".join(sorted(candidates))})
    _write_tsv(ambiguity_path, ["identifier_type", "identifier", "candidates"], rows)
    report = {
        "status": "complete",
        **mapper.qc(),
        "outputs": {"ambiguous_identifiers": _file_record(ambiguity_path)},
    }
    _write_json(output_dir / "hgnc_mapper_qc.json", report)
    return report


def _mapping_row(decision: GeneMappingDecision, **extra: Any) -> Dict[str, Any]:
    # Call-site annotations (for example duplicate-canonical rejection) must
    # be able to override the mapper's base decision without mutating it.
    return {**decision.row(), **extra}


def build_string_ppi(
    mapper: HGNCMapper,
    links_path: Path,
    info_path: Path,
    output_dir: Path,
    score_threshold: int = 700,
) -> Dict[str, Any]:
    if not 0 <= score_threshold <= 1000:
        raise ValueError("STRING score threshold must be between 0 and 1000")
    links_path = _require_file(Path(links_path), "STRING links")
    info_path = _require_file(Path(info_path), "STRING protein info")
    inputs = {"links": _file_record(links_path), "protein_info": _file_record(info_path)}
    output_dir.mkdir(parents=True, exist_ok=True)

    string_to_gene: Dict[str, str] = {}
    seen_string_ids: Set[str] = set()
    protein_rows: List[Dict[str, Any]] = []
    protein_status = Counter()
    with _open_text(info_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"#string_protein_id", "preferred_name"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"STRING protein.info missing columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            string_id = str(row.get("#string_protein_id", "")).strip()
            preferred = str(row.get("preferred_name", "")).strip()
            if not string_id:
                decision = GeneMappingDecision(
                    preferred, "", "", "", "rejected", "missing_string_id", (),
                    "STRING protein.info row has no identifier",
                )
            elif string_id in seen_string_ids:
                decision = GeneMappingDecision(
                    preferred, "", "", "", "rejected", "duplicate_string_id", (),
                    "duplicate STRING identifier rejected; first occurrence retained",
                )
            elif not string_id.startswith("9606."):
                decision = GeneMappingDecision(
                    preferred, "", "", "", "rejected", "non_human_string_id", (),
                    "STRING identifier does not have human taxid 9606",
                )
            else:
                decision = mapper.map(preferred)
            if string_id:
                seen_string_ids.add(string_id)
            protein_status[decision.method] += 1
            if decision.status == "mapped":
                string_to_gene[string_id] = decision.canonical_symbol
            protein_rows.append(_mapping_row(
                decision, row_number=row_number, string_protein_id=string_id, preferred_name=preferred,
            ))

    protein_map_path = output_dir / "string_protein_map.tsv.gz"
    _write_tsv(
        protein_map_path,
        ["row_number", "string_protein_id", "preferred_name", "source_symbol", "source_entrez_id",
         "canonical_symbol", "hgnc_id", "status", "method", "candidates", "reason"],
        protein_rows,
    )

    graph = nx.Graph()
    link_counts = Counter()
    rejection_path = output_dir / "string_edge_rejections.tsv.gz"
    with _atomic_gzip_text(rejection_path) as rejection_handle:
        rejection_writer = csv.DictWriter(
            rejection_handle,
            fieldnames=["line_number", "protein1", "protein2", "combined_score", "reason"],
            delimiter="\t",
        )
        rejection_writer.writeheader()
        with _open_text(links_path) as handle:
            header_line = handle.readline()
            if not header_line:
                raise ValueError("STRING links file is empty")
            header = header_line.strip().split()
            required = {"protein1", "protein2", "combined_score"}
            missing = required - set(header)
            if missing:
                raise ValueError(f"STRING links missing columns: {sorted(missing)}")
            p1_index, p2_index, score_index = (
                header.index("protein1"), header.index("protein2"), header.index("combined_score")
            )
            needed_length = max(p1_index, p2_index, score_index) + 1
            for line_number, line in enumerate(handle, start=2):
                link_counts["input_edges"] += 1
                parts = line.strip().split()
                if len(parts) < needed_length:
                    link_counts["malformed"] += 1
                    rejection_writer.writerow({"line_number": line_number, "reason": "malformed_row"})
                    continue
                protein1, protein2, raw_score = parts[p1_index], parts[p2_index], parts[score_index]
                try:
                    score = int(raw_score)
                except ValueError:
                    link_counts["invalid_score"] += 1
                    rejection_writer.writerow({
                        "line_number": line_number, "protein1": protein1, "protein2": protein2,
                        "combined_score": raw_score, "reason": "invalid_combined_score",
                    })
                    continue
                if score < 0 or score > 1000:
                    link_counts["out_of_range_score"] += 1
                    rejection_writer.writerow({
                        "line_number": line_number, "protein1": protein1, "protein2": protein2,
                        "combined_score": score, "reason": "combined_score_outside_0_to_1000",
                    })
                    continue
                if score < score_threshold:
                    link_counts["below_threshold"] += 1
                    continue
                link_counts["at_or_above_threshold"] += 1
                if not protein1.startswith("9606.") or not protein2.startswith("9606."):
                    link_counts["non_human"] += 1
                    rejection_writer.writerow({
                        "line_number": line_number, "protein1": protein1, "protein2": protein2,
                        "combined_score": score, "reason": "non_human_endpoint",
                    })
                    continue
                gene1, gene2 = string_to_gene.get(protein1, ""), string_to_gene.get(protein2, "")
                if not gene1 or not gene2:
                    reason = "both_endpoints_unmapped" if not gene1 and not gene2 else (
                        "protein1_unmapped" if not gene1 else "protein2_unmapped"
                    )
                    link_counts[reason] += 1
                    rejection_writer.writerow({
                        "line_number": line_number, "protein1": protein1, "protein2": protein2,
                        "combined_score": score, "reason": reason,
                    })
                    continue
                if gene1 == gene2:
                    link_counts["self_loop_after_mapping"] += 1
                    rejection_writer.writerow({
                        "line_number": line_number, "protein1": protein1, "protein2": protein2,
                        "combined_score": score, "reason": "self_loop_after_hgnc_mapping",
                    })
                    continue
                confidence = score / 1000.0
                if graph.has_edge(gene1, gene2):
                    link_counts["duplicate_after_mapping"] += 1
                    if score > int(graph[gene1][gene2]["combined_score"]):
                        graph[gene1][gene2].update(
                            combined_score=score, confidence=confidence, weight=confidence,
                        )
                        link_counts["duplicate_replaced_by_higher_score"] += 1
                    continue
                for string_id, canonical in ((protein1, gene1), (protein2, gene2)):
                    if canonical not in graph:
                        graph.add_node(
                            canonical,
                            hgnc_id=mapper.hgnc_ids.get(canonical, ""),
                            source="STRING",
                            source_string_id=string_id,
                        )
                graph.add_edge(
                    gene1, gene2, combined_score=score, confidence=confidence, weight=confidence,
                )
                link_counts["retained_unique_edges"] += 1

    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()
    largest = max((len(component) for component in nx.connected_components(graph)), default=0)
    graph.graph["metadata"] = {
        "source": "STRING",
        "version": STRING_RELEASE,
        "taxid": "9606",
        "identifier_namespace": "hgnc_symbol",
        "mapping_version": f"HGNC-{HGNC_RELEASE}",
        "score_threshold": score_threshold,
        "source_links_sha256": inputs["links"]["sha256"],
        "source_protein_info_sha256": inputs["protein_info"]["sha256"],
        "hgnc_sha256": mapper.sha256,
    }
    output_path = output_dir / PPI_OUTPUT_NAME
    deterministic_edges = []
    for left, right, attributes in graph.edges(data=True):
        gene_a, gene_b = sorted((str(left), str(right)))
        deterministic_edges.append({
            "gene_a": gene_a,
            "gene_b": gene_b,
            "combined_score": int(attributes["combined_score"]),
        })
    deterministic_edges.sort(key=lambda row: (
        row["gene_a"], row["gene_b"], -int(row["combined_score"]),
    ))
    _write_tsv(
        output_path,
        ["gene_a", "gene_b", "combined_score"],
        deterministic_edges,
    )

    for record in inputs.values():
        _assert_unchanged(record)
    report = {
        "status": "complete",
        "source": "STRING",
        "release": STRING_RELEASE,
        "score_threshold": score_threshold,
        "mapping_version": f"HGNC-{HGNC_RELEASE}",
        "inputs": inputs,
        "protein_mapping_counts": dict(sorted(protein_status.items())),
        "edge_counts": dict(sorted(link_counts.items())),
        "network": {
            "nodes": node_count,
            "edges": edge_count,
            "largest_connected_component_nodes": largest,
            "largest_connected_component_fraction": largest / max(1, node_count),
        },
        "outputs": {
            "network": _file_record(output_path),
            "safe_edge_table": _file_record(output_path),
            "protein_map": _file_record(protein_map_path),
            "edge_rejections": _file_record(rejection_path),
        },
        "notes": [
            "below-threshold rows are counted as an explicit filter, not copied to the rejection sidecar",
            "all malformed, non-human, unmapped-endpoint, and post-mapping self-loop rows are written",
            "duplicate mapped edges retain the maximum combined_score",
            "the primary PPI artifact is a deterministic gzip TSV edge table and never requires executable deserialization",
        ],
    }
    _write_json(output_dir / "string_ppi_qc.json", report)
    return report


def _normalize_cell_token(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value or "").strip().upper())


def _is_yes(value: Any) -> bool:
    return str(value or "").strip().lower() in {"yes", "true", "1", "y"}


def _read_score_cell_lines(combo_path: Path, expected_count: Optional[int] = 81) -> List[str]:
    delimiter = "\t" if combo_path.suffix.lower() in {".tsv", ".txt"} else ","
    seen: Set[str] = set()
    values: List[str] = []
    with combo_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if "cell_line_name" not in (reader.fieldnames or []):
            raise ValueError("ScoreRI combinations file requires cell_line_name")
        for row in reader:
            value = str(row.get("cell_line_name", "")).strip()
            if value and value not in seen:
                seen.add(value)
                values.append(value)
    values.sort(key=lambda value: (value.upper(), value))
    if expected_count is not None and len(values) != expected_count:
        raise ValueError(f"expected {expected_count} ScoreRI cell lines, found {len(values)}")
    return values


def _read_models(model_path: Path) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Set[str]]]:
    by_id: Dict[str, Dict[str, str]] = {}
    token_index: Dict[str, Set[str]] = defaultdict(set)
    with model_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"ModelID", "CellLineName", "StrippedCellLineName", "CCLEName"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"DepMap Model.csv missing columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            model_id = str(row.get("ModelID", "")).strip()
            if not model_id:
                continue
            if model_id in by_id:
                raise ValueError(f"duplicate ModelID in Model.csv at row {row_number}: {model_id}")
            cleaned = {key: str(value or "").strip() for key, value in row.items()}
            by_id[model_id] = cleaned
            candidates = [cleaned.get("CellLineName", ""), cleaned.get("StrippedCellLineName", "")]
            ccle = cleaned.get("CCLEName", "")
            if ccle:
                # Use the full CCLE identity only.  Taking the prefix before
                # ``_`` would incorrectly make resistant derivatives such as
                # A375_*_RESISTANT collide with the parental A-375 model.
                candidates.append(ccle)
            for candidate in candidates:
                token = _normalize_cell_token(candidate)
                if token:
                    token_index[token].add(model_id)
    return by_id, token_index


def _map_score_cell_lines(
    score_cell_lines: Sequence[str],
    models: Mapping[str, Mapping[str, str]],
    token_index: Mapping[str, Set[str]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for score_name in score_cell_lines:
        lookup_key = score_name.strip().upper()
        base: Dict[str, Any] = {
            "score_cell_line": score_name,
            "normalized_score_name": _normalize_cell_token(score_name),
            "model_id": "",
            "depmap_cell_line_name": "",
            "depmap_stripped_name": "",
            "mapping_status": "",
            "mapping_method": "",
            "candidate_model_ids": "",
            "reason": "",
            "expression_status": "not_scanned",
            "sequencing_id": "",
        }
        if lookup_key in EXPLICIT_SCORE_CELL_LINE_REJECTIONS:
            base.update({
                "mapping_status": "rejected",
                "mapping_method": "explicit_do_not_guess",
                "reason": EXPLICIT_SCORE_CELL_LINE_REJECTIONS[lookup_key],
            })
            rows.append(base)
            continue
        explicit = EXPLICIT_SCORE_CELL_LINE_ALIASES.get(lookup_key)
        if explicit:
            model_id = explicit["model_id"]
            model = models.get(model_id)
            if not model:
                base.update({
                    "mapping_status": "rejected", "mapping_method": "explicit_alias_missing_model",
                    "candidate_model_ids": model_id,
                    "reason": f"explicit model ID is absent from Model.csv: {model_id}",
                })
            elif str(model.get("CellLineName", "")).strip() != explicit["expected_depmap_name"]:
                base.update({
                    "mapping_status": "rejected", "mapping_method": "explicit_alias_identity_changed",
                    "candidate_model_ids": model_id,
                    "reason": (
                        f"explicit model identity mismatch: expected {explicit['expected_depmap_name']}, "
                        f"observed {model.get('CellLineName', '')}"
                    ),
                })
            else:
                base.update({
                    "model_id": model_id,
                    "depmap_cell_line_name": model.get("CellLineName", ""),
                    "depmap_stripped_name": model.get("StrippedCellLineName", ""),
                    "mapping_status": "mapped",
                    "mapping_method": "explicit_reviewed_alias",
                    "candidate_model_ids": model_id,
                    "reason": explicit["reason"],
                })
            rows.append(base)
            continue
        candidates = sorted(token_index.get(_normalize_cell_token(score_name), set()))
        base["candidate_model_ids"] = "|".join(candidates)
        if len(candidates) == 1:
            model_id = candidates[0]
            model = models[model_id]
            base.update({
                "model_id": model_id,
                "depmap_cell_line_name": model.get("CellLineName", ""),
                "depmap_stripped_name": model.get("StrippedCellLineName", ""),
                "mapping_status": "mapped",
                "mapping_method": "unique_normalized_exact",
                "reason": "unique match after case/punctuation normalization",
            })
        elif len(candidates) > 1:
            base.update({
                "mapping_status": "rejected", "mapping_method": "ambiguous_normalized_exact",
                "reason": "normalized name matches multiple DepMap ModelIDs",
            })
        else:
            base.update({
                "mapping_status": "rejected", "mapping_method": "unresolved",
                "reason": "no exact or explicitly reviewed DepMap model mapping",
            })
        rows.append(base)
    return rows


def _default_rna_profiles(profile_path: Optional[Path]) -> Dict[str, Set[str]]:
    if profile_path is None:
        return {}
    result: Dict[str, Set[str]] = defaultdict(set)
    with profile_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"ModelID", "SequencingID", "DataType", "IsDefaultEntryForModel"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"OmicsProfiles.csv missing columns: {sorted(missing)}")
        for row in reader:
            if str(row.get("DataType", "")).strip().lower() != "rna":
                continue
            if not _is_yes(row.get("IsDefaultEntryForModel")):
                continue
            model_id = str(row.get("ModelID", "")).strip()
            sequencing_id = str(row.get("SequencingID", "")).strip()
            if model_id and sequencing_id:
                result[model_id].add(sequencing_id)
    return result


_DEPMAP_GENE_HEADER = re.compile(r"^(.*?)\s+\((\d+)\)$")


def _parse_depmap_gene_header(value: str) -> Tuple[str, str]:
    match = _DEPMAP_GENE_HEADER.fullmatch(str(value).strip())
    if not match:
        return str(value).strip(), ""
    return match.group(1).strip(), match.group(2)


def build_depmap_expression(
    mapper: HGNCMapper,
    expression_path: Path,
    model_path: Path,
    score_combo_path: Path,
    profile_path: Optional[Path],
    output_dir: Path,
    expected_score_cell_lines: Optional[int] = 81,
) -> Dict[str, Any]:
    expression_path = _require_file(Path(expression_path), "DepMap expression")
    model_path = _require_file(Path(model_path), "DepMap Model.csv")
    score_combo_path = _require_file(Path(score_combo_path), "ScoreRI combinations")
    profile_path = _require_file(Path(profile_path), "DepMap OmicsProfiles.csv") if profile_path else None
    inputs = {
        "expression": _file_record(expression_path),
        "model": _file_record(model_path),
        "score_combinations": _file_record(score_combo_path),
    }
    if profile_path:
        inputs["omics_profiles"] = _file_record(profile_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    score_cell_lines = _read_score_cell_lines(score_combo_path, expected_score_cell_lines)
    models, token_index = _read_models(model_path)
    mapping_rows = _map_score_cell_lines(score_cell_lines, models, token_index)
    requested_model_ids = {
        str(row["model_id"]) for row in mapping_rows
        if row["mapping_status"] == "mapped" and row["model_id"]
    }
    profile_defaults = _default_rna_profiles(profile_path)

    with expression_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("DepMap expression file is empty") from exc
        required_meta = {"SequencingID", "ModelID", "IsDefaultEntryForModel"}
        missing = required_meta - set(header)
        if missing:
            raise ValueError(f"DepMap expression missing columns: {sorted(missing)}")
        sequencing_index = header.index("SequencingID")
        model_index = header.index("ModelID")
        default_index = header.index("IsDefaultEntryForModel")
        metadata_names = {
            "", "SequencingID", "ModelConditionID", "ModelID",
            "IsDefaultEntryForMC", "IsDefaultEntryForModel",
        }
        gene_indices = [index for index, name in enumerate(header) if name not in metadata_names]
        gene_rows: List[Dict[str, Any]] = []
        decisions: Dict[int, GeneMappingDecision] = {}
        canonical_to_indices: Dict[str, List[int]] = defaultdict(list)
        for index in gene_indices:
            source_symbol, entrez = _parse_depmap_gene_header(header[index])
            decision = mapper.map(source_symbol, entrez)
            decisions[index] = decision
            if decision.status == "mapped":
                canonical_to_indices[decision.canonical_symbol].append(index)
        duplicate_rejections: Dict[int, Tuple[str, str]] = {}
        included_mapped_indices: Set[int] = set()

        def mapping_priority(decision: GeneMappingDecision) -> int:
            if decision.method.startswith("entrez_unique"):
                return 0
            if decision.method == "approved_exact":
                return 1
            return 2

        for canonical, indices in canonical_to_indices.items():
            if len(indices) == 1:
                included_mapped_indices.add(indices[0])
                continue
            best_priority = min(mapping_priority(decisions[index]) for index in indices)
            best = [index for index in indices if mapping_priority(decisions[index]) == best_priority]
            if len(best) == 1:
                included_mapped_indices.add(best[0])
                for index in indices:
                    if index != best[0]:
                        duplicate_rejections[index] = (
                            "duplicate_canonical_lower_priority",
                            f"duplicate canonical symbol {canonical}; rejected in favor of unique higher-priority mapping",
                        )
            else:
                for index in indices:
                    duplicate_rejections[index] = (
                        "duplicate_canonical_equal_priority",
                        f"duplicate canonical symbol {canonical} has no unique highest-priority source column",
                    )
        included_gene_indices: List[int] = []
        for index in gene_indices:
            decision = decisions[index]
            included = decision.status == "mapped" and index in included_mapped_indices
            reason = decision.reason
            status = decision.status
            method = decision.method
            if index in duplicate_rejections:
                included = False
                status = "rejected"
                method, reason = duplicate_rejections[index]
            if included:
                included_gene_indices.append(index)
            gene_rows.append(_mapping_row(
                decision,
                source_column=header[index],
                source_column_index=index,
                status=status,
                method=method,
                included="yes" if included else "no",
                reason=reason,
            ))

        selected_rows: Dict[str, List[Tuple[int, str, List[str]]]] = defaultdict(list)
        scan_counts = Counter()
        for row_number, row in enumerate(reader, start=2):
            scan_counts["matrix_rows"] += 1
            if len(row) != len(header):
                scan_counts["malformed_matrix_rows"] += 1
                continue
            model_id = row[model_index].strip()
            if model_id not in requested_model_ids:
                continue
            scan_counts["requested_model_rows"] += 1
            if not _is_yes(row[default_index]):
                scan_counts["requested_nondefault_rows_filtered"] += 1
                continue
            sequencing_id = row[sequencing_index].strip()
            if profile_path and sequencing_id not in profile_defaults.get(model_id, set()):
                scan_counts["default_flag_crosscheck_rejected"] += 1
                continue
            selected_rows[model_id].append((row_number, sequencing_id, row))

    usable_profiles: Dict[str, Tuple[str, List[str]]] = {}
    for mapping in mapping_rows:
        if mapping["mapping_status"] != "mapped":
            mapping["expression_status"] = "unavailable_unresolved_cell_line"
            continue
        model_id = str(mapping["model_id"])
        profiles = selected_rows.get(model_id, [])
        if len(profiles) == 1:
            _, sequencing_id, values = profiles[0]
            usable_profiles[model_id] = (sequencing_id, values)
            mapping["expression_status"] = "available_default_profile"
            mapping["sequencing_id"] = sequencing_id
        elif not profiles:
            mapping["expression_status"] = "unavailable_no_verified_default_profile"
            mapping["reason"] = f"{mapping['reason']}; no verified default expression profile".strip("; ")
        else:
            mapping["expression_status"] = "unavailable_ambiguous_default_profiles"
            mapping["reason"] = (
                f"{mapping['reason']}; multiple verified default expression profiles: "
                + "|".join(sorted(profile[1] for profile in profiles))
            ).strip("; ")

    selected_mappings = [
        row for row in mapping_rows
        if row["expression_status"] == "available_default_profile"
    ]
    if not selected_mappings:
        raise ValueError("no ScoreRI cell line has a verified default DepMap expression profile")

    gene_map_path = output_dir / "gene_map.tsv"
    _write_tsv(
        gene_map_path,
        ["source_column_index", "source_column", "source_symbol", "source_entrez_id",
         "canonical_symbol", "hgnc_id", "status", "method", "candidates", "included", "reason"],
        gene_rows,
    )
    cell_map_path = output_dir / "cell_line_model_map.tsv"
    cell_fields = [
        "score_cell_line", "normalized_score_name", "model_id", "depmap_cell_line_name",
        "depmap_stripped_name", "mapping_status", "mapping_method", "candidate_model_ids",
        "expression_status", "sequencing_id", "reason",
    ]
    _write_tsv(cell_map_path, cell_fields, mapping_rows)
    cell_rejection_path = output_dir / "cell_line_mapping_rejections.tsv"
    _write_tsv(
        cell_rejection_path,
        cell_fields,
        [row for row in mapping_rows if row["expression_status"] != "available_default_profile"],
    )

    output_path = output_dir / EXPRESSION_OUTPUT_NAME
    missing_values = Counter()
    with _atomic_gzip_text(output_path) as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene"] + [row["score_cell_line"] for row in selected_mappings])
        for gene_index in included_gene_indices:
            canonical = decisions[gene_index].canonical_symbol
            output_row: List[str] = [canonical]
            for mapping in selected_mappings:
                model_id = str(mapping["model_id"])
                raw_value = usable_profiles[model_id][1][gene_index].strip()
                try:
                    value = float(raw_value)
                except ValueError:
                    value = float("nan")
                if not math.isfinite(value):
                    output_row.append("")
                    missing_values[str(mapping["score_cell_line"])] += 1
                else:
                    # The source is already log2(TPM+1); preserve its value and
                    # never apply another transformation.
                    output_row.append(raw_value)
            writer.writerow(output_row)

    for record in inputs.values():
        _assert_unchanged(record)
    mapping_counts = Counter(row["mapping_method"] for row in mapping_rows)
    expression_counts = Counter(row["expression_status"] for row in mapping_rows)
    gene_counts = Counter(row["method"] for row in gene_rows)
    report = {
        "status": "complete" if len(selected_mappings) == len(score_cell_lines) else "complete_with_declared_gaps",
        "source": "DepMap",
        "release": DEPMAP_RELEASE,
        "value_semantics": "source values preserved as log2(TPM+1); no additional log transformation",
        "mapping_version": f"HGNC-{HGNC_RELEASE}",
        "inputs": inputs,
        "score_cell_lines": {
            "requested": len(score_cell_lines),
            "with_verified_default_expression": len(selected_mappings),
            "mapping_method_counts": dict(sorted(mapping_counts.items())),
            "expression_status_counts": dict(sorted(expression_counts.items())),
            "explicit_unresolved": [
                row["score_cell_line"] for row in mapping_rows
                if row["mapping_method"] == "explicit_do_not_guess"
            ],
        },
        "genes": {
            "source_columns": len(gene_indices),
            "included_canonical_genes": len(included_gene_indices),
            "duplicate_canonical_columns_rejected": len(duplicate_rejections),
            "duplicate_canonical_groups_resolved_by_mapping_priority": sum(
                len(indices) > 1 and any(index in included_mapped_indices for index in indices)
                for indices in canonical_to_indices.values()
            ),
            "mapping_method_counts": dict(sorted(gene_counts.items())),
        },
        "matrix_scan_counts": dict(sorted(scan_counts.items())),
        "missing_or_nonfinite_values_by_cell_line": dict(sorted(missing_values.items())),
        "outputs": {
            "expression": _file_record(output_path),
            "cell_line_model_map": _file_record(cell_map_path),
            "cell_line_rejections": _file_record(cell_rejection_path),
            "gene_map": _file_record(gene_map_path),
        },
    }
    _write_json(output_dir / "depmap_expression_qc.json", report)
    return report


def build_reactome(
    mapper: HGNCMapper,
    reactome_path: Path,
    output_dir: Path,
) -> Dict[str, Any]:
    reactome_path = _require_file(Path(reactome_path), "Reactome GMT")
    inputs = {"reactome_gmt": _file_record(reactome_path)}
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / REACTOME_OUTPUT_NAME
    rejection_path = output_dir / "reactome_gene_rejections.tsv.gz"
    pathway_rejection_path = output_dir / "reactome_pathway_rejections.tsv"
    counts = Counter()
    rejection_rows: List[Dict[str, Any]] = []
    pathway_rejections: List[Dict[str, Any]] = []
    with _atomic_text(output_path, newline="\n") as output_handle:
        with reactome_path.open("r", encoding="utf-8", errors="strict") as input_handle:
            for line_number, line in enumerate(input_handle, start=1):
                counts["input_pathways"] += 1
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 3:
                    counts["malformed_pathways"] += 1
                    pathway_rejections.append({
                        "line_number": line_number,
                        "pathway_name": parts[0] if parts else "",
                        "reason": "GMT row has fewer than three fields",
                    })
                    continue
                pathway_name, stable_id = parts[0].strip(), parts[1].strip()
                if not pathway_name or not stable_id:
                    counts["malformed_pathways"] += 1
                    pathway_rejections.append({
                        "line_number": line_number, "pathway_name": pathway_name,
                        "reason": "pathway name or Reactome stable ID is blank",
                    })
                    continue
                canonical_members: List[str] = []
                seen: Set[str] = set()
                for source_gene in parts[2:]:
                    counts["input_memberships"] += 1
                    decision = mapper.map(source_gene)
                    if decision.status != "mapped":
                        counts[f"rejected_{decision.method}"] += 1
                        rejection_rows.append(_mapping_row(
                            decision,
                            line_number=line_number,
                            pathway_name=pathway_name,
                            reactome_stable_id=stable_id,
                        ))
                        continue
                    if decision.method != "approved_exact":
                        counts["canonicalized_alias_memberships"] += 1
                    if decision.canonical_symbol in seen:
                        counts["duplicate_memberships_after_mapping"] += 1
                        continue
                    seen.add(decision.canonical_symbol)
                    canonical_members.append(decision.canonical_symbol)
                if not canonical_members:
                    counts["empty_pathways_after_mapping"] += 1
                    pathway_rejections.append({
                        "line_number": line_number, "pathway_name": pathway_name,
                        "reason": "no HGNC-approved human members remain",
                    })
                    continue
                output_handle.write("\t".join([pathway_name, stable_id] + canonical_members) + "\n")
                counts["retained_pathways"] += 1
                counts["retained_memberships"] += len(canonical_members)

    _write_tsv(
        rejection_path,
        ["line_number", "pathway_name", "reactome_stable_id", "source_symbol", "source_entrez_id",
         "canonical_symbol", "hgnc_id", "status", "method", "candidates", "reason"],
        rejection_rows,
    )
    _write_tsv(pathway_rejection_path, ["line_number", "pathway_name", "reason"], pathway_rejections)
    _assert_unchanged(inputs["reactome_gmt"])
    report = {
        "status": "complete",
        "source": "Reactome",
        "release": REACTOME_RELEASE,
        "taxid": "9606",
        "mapping_version": f"HGNC-{HGNC_RELEASE}",
        "inputs": inputs,
        "counts": dict(sorted(counts.items())),
        "outputs": {
            "cleaned_gmt": _file_record(output_path),
            "gene_rejections": _file_record(rejection_path),
            "pathway_rejections": _file_record(pathway_rejection_path),
        },
    }
    _write_json(output_dir / "reactome_qc.json", report)
    return report


def build_cosmic(
    mapper: HGNCMapper,
    cosmic_path: Path,
    output_dir: Path,
) -> Dict[str, Any]:
    cosmic_path = _require_file(Path(cosmic_path), "COSMIC Cancer Gene Census")
    inputs = {"cosmic_cgc": _file_record(cosmic_path)}
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / COSMIC_OUTPUT_NAME
    tier1_path = output_dir / COSMIC_TIER1_OUTPUT_NAME
    rejection_path = output_dir / "cosmic_gene_rejections.tsv"
    counts = Counter()
    rejections: List[Dict[str, Any]] = []
    tier1: Set[str] = set()

    with cosmic_path.open("r", encoding="utf-8-sig", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        required = {"GENE_SYMBOL", "TIER", "ROLE_IN_CANCER", "TUMOUR_TYPES_SOMATIC"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"COSMIC CGC missing columns: {sorted(missing)}")
        original_fields = list(reader.fieldnames or [])
        audit_fields = ["SOURCE_GENE_SYMBOL", "HGNC_ID", "HGNC_MAPPING_METHOD"]
        with _atomic_text(output_path) as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=original_fields + audit_fields,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            for row_number, row in enumerate(reader, start=2):
                counts["input_rows"] += 1
                source_symbol = str(row.get("GENE_SYMBOL", "")).strip()
                decision = mapper.map(source_symbol)
                if decision.status != "mapped":
                    counts[f"rejected_{decision.method}"] += 1
                    rejections.append(_mapping_row(decision, row_number=row_number))
                    continue
                if decision.method != "approved_exact":
                    counts["canonicalized_alias_rows"] += 1
                output_row = dict(row)
                output_row["SOURCE_GENE_SYMBOL"] = source_symbol
                output_row["GENE_SYMBOL"] = decision.canonical_symbol
                output_row["HGNC_ID"] = decision.hgnc_id
                output_row["HGNC_MAPPING_METHOD"] = decision.method
                writer.writerow(output_row)
                counts["retained_rows"] += 1
                if str(row.get("TIER", "")).strip() == "1":
                    tier1.add(decision.canonical_symbol)
                elif str(row.get("TIER", "")).strip() not in {"", "2"}:
                    counts["unexpected_tier_values"] += 1

    if not tier1:
        raise ValueError("COSMIC CGC yielded no canonical pan-cancer tier-1 genes")
    with _atomic_text(tier1_path, newline="\n") as handle:
        for symbol in sorted(tier1):
            handle.write(symbol + "\n")
    _write_tsv(
        rejection_path,
        ["row_number", "source_symbol", "source_entrez_id", "canonical_symbol", "hgnc_id",
         "status", "method", "candidates", "reason"],
        rejections,
    )
    _assert_unchanged(inputs["cosmic_cgc"])
    report = {
        "status": "complete",
        "source": "COSMIC Cancer Gene Census",
        "release": COSMIC_RELEASE,
        "genome_build": "GRCh38",
        "taxid": "9606",
        "mapping_version": f"HGNC-{HGNC_RELEASE}",
        "inputs": inputs,
        "counts": {**dict(sorted(counts.items())), "pan_cancer_tier1_genes": len(tier1)},
        "outputs": {
            "canonical_tsv": _file_record(output_path),
            "pan_cancer_tier1_gene_set": _file_record(tier1_path),
            "gene_rejections": _file_record(rejection_path),
        },
    }
    _write_json(output_dir / "cosmic_qc.json", report)
    return report


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build HGNC-normalized production context artifacts without network access"
    )
    parser.add_argument("command", choices=("hgnc", "ppi", "depmap", "reactome", "cosmic", "all"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--hgnc", type=Path, default=DEFAULT_HGNC)
    parser.add_argument("--string-links", type=Path, default=DEFAULT_STRING_LINKS)
    parser.add_argument("--string-info", type=Path, default=DEFAULT_STRING_INFO)
    parser.add_argument("--string-score-threshold", type=int, default=700)
    parser.add_argument("--depmap-expression", type=Path, default=DEFAULT_DEPMAP_EXPRESSION)
    parser.add_argument("--depmap-model", type=Path, default=DEFAULT_DEPMAP_MODEL)
    parser.add_argument("--depmap-profiles", type=Path, default=DEFAULT_DEPMAP_PROFILES)
    parser.add_argument("--score-combos", type=Path, default=DEFAULT_SCORE_COMBOS)
    parser.add_argument("--expected-score-cell-lines", type=int, default=81)
    parser.add_argument("--reactome", type=Path, default=DEFAULT_REACTOME)
    parser.add_argument("--cosmic", type=Path, default=DEFAULT_COSMIC)
    return parser


def run_build(args: argparse.Namespace, verify_hgnc_checksum: bool = True) -> Dict[str, Any]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    mapper = HGNCMapper(
        Path(args.hgnc),
        expected_sha256=HGNC_EXPECTED_SHA256 if verify_hgnc_checksum else None,
    )
    steps: Dict[str, Any] = {"hgnc": build_hgnc_reports(mapper, output_dir)}
    commands = {args.command} if args.command != "all" else {"ppi", "depmap", "reactome", "cosmic"}
    if "ppi" in commands:
        steps["ppi"] = build_string_ppi(
            mapper, Path(args.string_links), Path(args.string_info), output_dir,
            score_threshold=int(args.string_score_threshold),
        )
    if "depmap" in commands:
        expected = int(args.expected_score_cell_lines)
        if expected <= 0:
            raise ValueError("--expected-score-cell-lines must be positive")
        steps["depmap"] = build_depmap_expression(
            mapper,
            Path(args.depmap_expression),
            Path(args.depmap_model),
            Path(args.score_combos),
            Path(args.depmap_profiles) if args.depmap_profiles else None,
            output_dir,
            expected_score_cell_lines=expected,
        )
    if "reactome" in commands:
        steps["reactome"] = build_reactome(mapper, Path(args.reactome), output_dir)
    if "cosmic" in commands:
        steps["cosmic"] = build_cosmic(mapper, Path(args.cosmic), output_dir)
    report = {
        "schema_version": "1.0",
        "builder": str(Path(__file__).resolve()),
        "command": args.command,
        "status": (
            "complete_with_declared_gaps"
            if any(step.get("status") == "complete_with_declared_gaps" for step in steps.values())
            else "complete"
        ),
        "completed_at_utc": _utc_now(),
        "fixed_releases": {
            "HGNC": HGNC_RELEASE,
            "STRING": STRING_RELEASE,
            "DepMap": DEPMAP_RELEASE,
            "Reactome": REACTOME_RELEASE,
            "COSMIC_CGC": COSMIC_RELEASE,
        },
        "steps": steps,
    }
    _write_json(output_dir / "build_qc.json", report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        report = run_build(args, verify_hgnc_checksum=True)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({
        "status": report["status"],
        "command": report["command"],
        "output_dir": str(Path(args.output_dir).resolve()),
        "build_qc": str((Path(args.output_dir) / "build_qc.json").resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
