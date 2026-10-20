"""Local + PLATO target resolution, caching, human filtering, and fusion.

PLATO remains an external prediction source.  This module consumes PLATO JSON
results, including audited one-request form responses, and keeps the raw
response bytes for provenance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from contextlib import closing
import csv
import gzip
import json
import math
import os
import re
import sqlite3
import tempfile


HUMAN_LABELS = {"homo sapiens", "human", "h sapiens"}
NON_PROTEIN_TARGET_TYPES = {
    "cell line", "organism", "tissue", "subcellular", "unknown", "non protein", "non molecular",
    "non protein target", "nucleic acid", "dna", "rna", "small molecule", "macromolecule",
}
PROTEIN_TARGET_TYPES = {
    "protein",
    "single protein",
    "protein complex",
    "protein family",
    "protein protein interaction",
    "chimeric protein",
    "protein nucleic acid complex",
}
PLATO_PROVIDER_LABELS = {"plato"}
PLATO_TARGET_FISHING_METHODS = {"target fishing", "targetfishing"}
PLATO_TARGET_FISHING_SCORE_RANGE = (0.0, 13.0)
VERIFIED_QUERY_BINDING_STATUSES = {
    "payload_query_structure_verified",
    "plato_single_query_post_verified",
}

# Production-approved PLATO schema fixtures.  This allowlist is deliberately
# empty: the tracked Gemcitabine r37 response remains a parser fixture
# candidate whose manifest says ``not_yet_allowlisted_for_scoring``.  Its
# observed field shape may be parsed for regression coverage, but it must not
# establish provider, release, method, or production schema trust until an
# independent fixture review is completed.
VERIFIED_PLATO_SCHEMA_FINGERPRINTS: Tuple[str, ...] = ()

PARSER_VERSION = "4.0"
CACHE_SCHEMA_VERSION = "3"


def _text(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def _first(record: Mapping[str, Any], names: Sequence[str], default: Any = "") -> Any:
    lowered = {str(key).strip().lower(): value for key, value in record.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value not in (None, ""):
            return value
    return default


def _bool(value: Any) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    normalized = _text(value).lower()
    if normalized in {"yes", "true", "1", "reliable", "high"}:
        return True
    if normalized in {"no", "not", "false", "0", "not reliable", "unreliable", "low"}:
        return False
    return None


def _float(value: Any) -> Optional[float]:
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError, OverflowError):
        return None


def _int(value: Any) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError, OverflowError):
        return None


def _label(value: Any) -> str:
    if isinstance(value, Mapping):
        return _text(_first(value, ("name", "label", "value", "text", "id")))
    return _text(value)


def _taxid(value: Any) -> str:
    label = _label(value)
    try:
        number = float(label)
        if math.isfinite(number) and number.is_integer():
            return str(int(number))
    except (TypeError, ValueError):
        pass
    prefixed = re.fullmatch(
        r"(?:ncbi\s*(?:taxonomy|taxon)|tax(?:onomy|on)?(?:\s*id)?)\s*[:#]?\s*(\d+)",
        label.strip(),
        flags=re.I,
    )
    if prefixed:
        return prefixed.group(1)
    return label


def _normalize_identifier(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", _text(value).upper())


def _normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", _text(value).lower()).strip()


def _looks_like_uniprot_accession(value: str) -> bool:
    """Recognize canonical UniProt accessions, including 10-character accessions."""
    normalized = _text(value).upper()
    return bool(re.fullmatch(
        r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?",
        normalized,
    ))


def _looks_like_gene_symbol(value: str) -> bool:
    """Conservatively recognize a symbol-like raw target identifier."""
    normalized = _text(value)
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9.-]{0,19}", normalized))


def _canonical_hgnc_symbol(value: Any) -> str:
    """Preserve the approved HGNC spelling used by biological context files.

    Identifier lookup keys may remove punctuation, but the symbol handed to
    PPI/expression/pathway code must retain meaningful characters such as the
    hyphen in ``HLA-A`` and ``MT-ND1``.
    """
    symbol = _text(value).upper()
    return symbol if _looks_like_gene_symbol(symbol) else ""


def _gene_symbol_lookup_key(value: Any) -> str:
    """Return a comparison-only key; never expose this as the HGNC symbol."""
    return _normalize_identifier(_text(value))


def _is_human_label(value: Any) -> bool:
    normalized = _normalize_name(_label(value))
    return normalized in HUMAN_LABELS


def _validated_plato_score(value: Any) -> Tuple[Optional[float], List[str]]:
    """Validate PLATO target-fishing's documented 0..13 ranking score."""
    if value in (None, ""):
        return None, ["plato_score_missing"]
    if isinstance(value, bool):
        return None, ["plato_score_invalid_type"]
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None, ["plato_score_not_numeric"]
    if not math.isfinite(number):
        return None, ["plato_score_non_finite"]
    lower, upper = PLATO_TARGET_FISHING_SCORE_RANGE
    if not lower <= number <= upper:
        return None, ["plato_score_out_of_range_0_13"]
    return number, []


def standardize_smiles(smiles: str) -> Dict[str, Any]:
    """Create an auditable structure identity; only RDKit identity is valid."""
    submitted = _text(smiles)
    result: Dict[str, Any] = {
        "submitted_smiles": submitted,
        "canonical_isomeric_smiles": "",
        "inchikey": "",
        "valid": False,
        "standardization_status": "not_validated",
        "warnings": [],
    }
    if not submitted:
        result["standardization_status"] = "invalid_empty"
        result["warnings"].append("empty_smiles")
        return result
    try:
        from rdkit import Chem

        molecule = Chem.MolFromSmiles(submitted)
        if molecule is None:
            result["standardization_status"] = "invalid_rdkit_parse"
            result["warnings"].append("rdkit_parse_failed")
            return result
        result["canonical_isomeric_smiles"] = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        try:
            result["inchikey"] = Chem.MolToInchiKey(molecule)
        except (AttributeError, RuntimeError, ValueError):
            result["warnings"].append("inchikey_unavailable")
        if result["canonical_isomeric_smiles"] and result["inchikey"]:
            result["valid"] = True
            result["standardization_status"] = "rdkit_canonical_identity"
        else:
            result["standardization_status"] = "invalid_incomplete_rdkit_identity"
    except ImportError:
        result["standardization_status"] = "validation_unavailable_rdkit_missing"
        result["warnings"].append("rdkit_unavailable")
    except Exception as exc:  # RDKit builds vary in their exception classes.
        result["standardization_status"] = "validation_error"
        result["warnings"].append(f"rdkit_standardization_error:{type(exc).__name__}")
    return result


def stereo_agnostic_smiles(smiles: str) -> str:
    """Remove only stereochemistry; retain isotopes, charge, bonds and components."""
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(_text(smiles))
    if molecule is None or not _text(smiles):
        raise ValueError("A valid nonempty SMILES is required for structure grouping")
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def structure_identity_matches(query: Mapping[str, Any], expected: Mapping[str, Any]) -> bool:
    """Compare identities under an explicit request policy without rewriting either."""
    actual = standardize_smiles(_text(query.get("canonical_isomeric_smiles")))
    if not actual.get("valid") or actual.get("inchikey") != query.get("inchikey"):
        return False
    if expected.get("ignore_stereochemistry") is True:
        try:
            return stereo_agnostic_smiles(actual["canonical_isomeric_smiles"]) == stereo_agnostic_smiles(
                _text(expected.get("canonical_isomeric_smiles")))
        except (ValueError, ImportError):
            return False
    return (actual["canonical_isomeric_smiles"] == expected.get("canonical_isomeric_smiles")
            and actual["inchikey"] == expected.get("inchikey"))


def _legacy_smiles_cache_key(smiles: str) -> str:
    return sha256(_text(smiles).encode("utf-8")).hexdigest()


def smiles_cache_key(smiles: str) -> str:
    """Hash canonical isomeric SMILES when RDKit is available."""
    structure = standardize_smiles(smiles)
    identity = structure["canonical_isomeric_smiles"] or structure["submitted_smiles"]
    return sha256(identity.encode("utf-8")).hexdigest()


def _payload_sha256(payload: Any) -> str:
    normalized = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return sha256(normalized.encode("utf-8")).hexdigest()


def _is_sha256(value: Any) -> bool:
    return bool(re.fullmatch(r"[0-9a-f]{64}", _text(value)))


def _payload_query_smiles(payload: Any) -> List[str]:
    """Read only query-level structure fields, never ligand/target records."""
    if not isinstance(payload, Mapping):
        return []
    candidates: List[str] = []
    direct_keys = (
        "query_smiles", "input_smiles", "submitted_smiles", "canonical_smiles",
        "canonical_isomeric_smiles", "smiles",
    )

    def collect(container: Any, *, allow_scalar: bool = False) -> None:
        if allow_scalar and isinstance(container, str) and _text(container):
            candidates.append(_text(container))
            return
        if not isinstance(container, Mapping):
            return
        for key in direct_keys:
            value = _first(container, (key,), "")
            if isinstance(value, str) and _text(value):
                candidates.append(_text(value))

    collect(payload)
    for key in ("metadata", "request", "input", "query", "molecule", "compound"):
        nested = _first(payload, (key,), None)
        collect(nested, allow_scalar=key == "query")
    return list(dict.fromkeys(candidates))


def _looks_like_plato_r37_target_fishing_payload(payload: Any) -> bool:
    """Recognize the raw JSON array returned by PLATO r37 target fishing."""
    if not isinstance(payload, list) or not payload:
        return False
    required = {"chemblid", "nome", "score", "simili", "tid", "validato"}
    records = [item for item in payload if isinstance(item, Mapping)]
    if len(records) != len(payload):
        return False
    for record in records:
        keys = {str(key).strip().lower() for key in record}
        if not required <= keys:
            return False
        if not isinstance(record.get("simili"), list):
            return False
    return True


def _query_binding_report(payload: Any, submitted_smiles: str) -> Dict[str, Any]:
    expected = standardize_smiles(submitted_smiles)
    candidates = _payload_query_smiles(payload)
    if not expected.get("valid"):
        return {
            "query_binding_status": "query_structure_not_validated",
            "payload_query_smiles": candidates,
            "query_binding_method": "unverified",
            "query_binding_warnings": list(expected.get("warnings", [])),
        }
    if not candidates:
        return {
            "query_binding_status": "user_asserted_unverified",
            "payload_query_smiles": [],
            "query_binding_method": "payload_query_structure_absent",
            "query_binding_warnings": ["payload_query_structure_absent"],
        }
    expected_identity = (expected["canonical_isomeric_smiles"], expected["inchikey"])
    mismatches: List[str] = []
    invalid_candidates: List[str] = []
    for candidate in candidates:
        observed = standardize_smiles(candidate)
        if not observed.get("valid"):
            invalid_candidates.append(candidate)
            continue
        observed_identity = (observed["canonical_isomeric_smiles"], observed["inchikey"])
        if observed_identity != expected_identity:
            mismatches.append(candidate)
    if mismatches:
        raise ValueError(
            "PLATO query structure mismatch: payload query SMILES does not match the submitted SMILES"
        )
    if invalid_candidates:
        return {
            "query_binding_status": "payload_query_structure_not_validated",
            "payload_query_smiles": candidates,
            "invalid_payload_query_smiles": invalid_candidates,
            "query_binding_method": "unverified",
            "query_binding_warnings": ["payload_query_structure_invalid_or_unavailable"],
        }
    return {
        "query_binding_status": "payload_query_structure_verified",
        "payload_query_smiles": candidates,
        "query_binding_method": "rdkit_canonical_smiles_and_inchikey",
        "query_binding_warnings": [],
    }


def _single_query_binding_report(
    submitted_smiles: str,
    raw_sha256: str,
    evidence: Optional[Mapping[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Verify a user-triggered single PLATO POST request against the cache raw bytes."""
    if not evidence:
        return None
    expected = standardize_smiles(submitted_smiles)
    if expected.get("valid") is not True:
        return None
    endpoint = _text(
        evidence.get("response_final_url") or evidence.get("endpoint")
    ).rstrip("/")
    requested_endpoint = _text(
        evidence.get("requested_endpoint") or evidence.get("endpoint")
    ).rstrip("/")
    method = _text(evidence.get("http_method") or evidence.get("method")).upper()
    response_sha256 = _text(evidence.get("response_sha256")).lower()
    form = evidence.get("form") if isinstance(evidence.get("form"), Mapping) else {}
    candidates = [
        _text(evidence.get("submitted_smiles")),
        _text(evidence.get("query_smiles")),
        _text(form.get("smiles")) if isinstance(form, Mapping) else "",
    ]
    candidates = [value for value in dict.fromkeys(candidates) if value]
    allowed_endpoint = endpoint == "https://prometheus.farmacia.uniba.it/plato/pred.html"
    allowed_requested_endpoint = (
        requested_endpoint == "https://prometheus.farmacia.uniba.it/plato/pred.html"
    )
    warnings: List[str] = []
    if not allowed_endpoint:
        warnings.append("single_query_endpoint_not_allowlisted")
    if not allowed_requested_endpoint:
        warnings.append("single_query_requested_endpoint_not_allowlisted")
    if method != "POST":
        warnings.append("single_query_method_not_post")
    if response_sha256 and response_sha256 != _text(raw_sha256).lower():
        warnings.append("single_query_response_sha256_mismatch")
    if not response_sha256:
        warnings.append("single_query_response_sha256_missing")
    if isinstance(form, Mapping):
        if _text(form.get("tipoout")).lower() != "json":
            warnings.append("single_query_output_not_json")
        if _text(form.get("tipopred")).lower() != "mussel":
            warnings.append("single_query_method_not_target_fishing")
    else:
        warnings.append("single_query_form_missing")
    if not candidates:
        warnings.append("single_query_smiles_missing")
    for candidate in candidates:
        observed = standardize_smiles(candidate)
        if (
            observed.get("valid") is not True
            or observed.get("canonical_isomeric_smiles") != expected.get("canonical_isomeric_smiles")
            or str(observed.get("inchikey", "")).upper() != str(expected.get("inchikey", "")).upper()
        ):
            warnings.append("single_query_smiles_identity_mismatch")
            break
    if warnings:
        return None
    return {
        "query_binding_status": "plato_single_query_post_verified",
        "payload_query_smiles": [],
        "query_binding_method": "audited_single_plato_form_post",
        "query_binding_warnings": [],
        "single_query_evidence": {
            "endpoint": endpoint,
            "requested_endpoint": requested_endpoint,
            "http_method": method,
            "submitted_smiles": candidates[0],
            "response_sha256": response_sha256,
            "form": {
                "tipoout": _text(form.get("tipoout")),
                "tipopred": _text(form.get("tipopred")),
            },
        },
    }


def _loads_strict_json(text: str) -> Any:
    return json.loads(
        text,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Non-finite JSON constant: {value}")),
    )


def _json_safe_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe_value(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (str, int, bool)):
        return value
    return str(value)


def _cache_mapping(
    value: Any,
    field_name: str,
    warnings: List[str],
) -> Dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        warnings.append(f"cached_parser_report_{field_name}_invalid_type")
        return {}
    return dict(value)


def _cache_string_list(
    value: Any,
    field_name: str,
    warnings: List[str],
) -> List[str]:
    if value is None:
        return []
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        warnings.append(f"cached_parser_report_{field_name}_invalid_type")
        return []
    return list(value)


def _cache_bool(value: Any, field_name: str, warnings: List[str]) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    if value is not None:
        warnings.append(f"cached_parser_report_{field_name}_invalid_type")
    return None


def _sqlite_resource_metadata(connection: sqlite3.Connection) -> Dict[str, str]:
    """Read optional artifact metadata without weakening older schema support."""
    try:
        exists = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='resource_metadata'"
        ).fetchone()
        if not exists:
            return {}
        columns = {
            str(row[1]).lower() for row in connection.execute("PRAGMA table_info(resource_metadata)")
        }
        if not {"key", "value"} <= columns:
            return {}
        return {
            _text(key): _text(value)
            for key, value in connection.execute("SELECT key, value FROM resource_metadata")
            if _text(key)
        }
    except sqlite3.Error:
        return {}


@dataclass
class PlatoTargetEvidence:
    target_name: str
    target_symbol: str = ""
    target_chembl_id: str = ""
    uniprot_id: str = ""
    organism: str = ""
    target_type: str = ""
    taxid: str = ""
    score: Optional[float] = None
    reliable: Optional[bool] = None
    reliable_raw: str = ""
    rank: Optional[int] = None
    applicability_score: Optional[float] = None
    action: str = ""
    score_semantics: str = "PLATO provider ranking score; not a calibrated probability"
    supporting_ligands: List[Dict[str, Any]] = field(default_factory=list)
    parse_warnings: List[str] = field(default_factory=list)
    raw: Dict[str, Any] = field(default_factory=dict)

    @property
    def target_id(self) -> str:
        return self.target_symbol or self.uniprot_id or self.target_chembl_id or self.target_name

    @property
    def alignment_keys(self) -> List[str]:
        keys: List[str] = []
        for prefix, value in (
            ("uniprot", self.uniprot_id),
            ("chembl", self.target_chembl_id),
            ("symbol", self.target_symbol),
        ):
            normalized = _normalize_identifier(value)
            if normalized:
                keys.append(f"{prefix}:{normalized}")
        name = _normalize_name(self.target_name)
        if name:
            keys.append(f"name:{name}")
        return keys

    @property
    def human_status(self) -> str:
        normalized_taxid = _taxid(self.taxid)
        normalized = _normalize_name(self.organism)
        if normalized_taxid:
            taxid_human = normalized_taxid == "9606"
            if normalized and taxid_human != _is_human_label(self.organism):
                return "species_conflict"
            return "human_explicit" if taxid_human else "non_human"
        if _is_human_label(self.organism):
            return "human_explicit"
        if normalized:
            return "non_human"
        return "unknown_species"

    def is_human(self) -> bool:
        return self.human_status.startswith("human_")

    def is_protein_target(self) -> bool:
        normalized = _normalize_name(self.target_type)
        if not normalized or normalized in NON_PROTEIN_TARGET_TYPES:
            return False
        return normalized in PROTEIN_TARGET_TYPES

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data.pop("raw", None)
        data["human_status"] = self.human_status
        return _json_safe_value(data)


@dataclass
class FusedTarget:
    target_id: str
    target_name: str
    source: str
    evidence_status: str
    local: Optional[Dict[str, Any]] = None
    plato: Optional[PlatoTargetEvidence] = None
    alignment_method: str = ""
    alignment_status: str = "unmatched"

    def to_tuple(self) -> Tuple[str, str]:
        return self.target_id, self.target_name

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target_id": self.target_id,
            "target_name": self.target_name,
            "source": self.source,
            "evidence_status": self.evidence_status,
            "local": _json_safe_value(self.local),
            "plato": self.plato.to_dict() if self.plato else None,
            "alignment_method": self.alignment_method,
            "alignment_status": self.alignment_status,
        }


class ChemblTargetMapper:
    """Translate PLATO ChEMBL target IDs into HGNC symbols for context mapping.

    This is not a target source.  It only enriches targets already predicted by
    PLATO so they can be matched to HGNC-normalized PPI/expression/pathway
    context.
    """

    REQUIRED_COLUMNS = {
        "target_chembl_id",
        "target_name",
        "chembl_target_type",
        "component_uniprot_id",
        "canonical_hgnc_symbol",
        "hgnc_id",
        "mapping_status",
        "target_component_count",
    }

    def __init__(
        self,
        path: str | Path | None = None,
        *,
        expected_sha256: str = "",
        source: str = "ChEMBL",
        release: str = "r37",
    ):
        self.path = Path(path) if path else None
        self.expected_sha256 = _text(expected_sha256).lower()
        self.source = source
        self.release = release
        self.mapping: Dict[str, Dict[str, Any]] = {}
        self.status = "disabled"
        self.sha256 = ""
        self.load_warnings: List[str] = []
        if self.path:
            self._load()

    @staticmethod
    def _file_sha256(path: Path) -> str:
        digest = sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _open_text(self):
        if self.path is None:
            raise FileNotFoundError("ChEMBL target mapping path is not configured")
        if self.path.suffix.lower() == ".gz":
            return gzip.open(self.path, "rt", encoding="utf-8", newline="")
        return self.path.open("r", encoding="utf-8-sig", newline="")

    def _load(self) -> None:
        if self.path is None:
            self.status = "disabled"
            return
        if not self.path.is_file():
            self.status = "missing"
            self.load_warnings.append("chembl_target_mapping_missing")
            return
        self.sha256 = self._file_sha256(self.path)
        if self.expected_sha256 and self.sha256.lower() != self.expected_sha256:
            self.status = "sha256_mismatch"
            self.load_warnings.append("chembl_target_mapping_sha256_mismatch")
            return
        grouped: Dict[str, List[Dict[str, Any]]] = {}
        invalid_symbols = 0
        with self._open_text() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = self.REQUIRED_COLUMNS - set(reader.fieldnames or [])
            if missing:
                self.status = "invalid_schema"
                self.load_warnings.append(f"chembl_target_mapping_missing_columns:{','.join(sorted(missing))}")
                return
            for row in reader:
                if _text(row.get("mapping_status")).lower() != "mapped":
                    continue
                target_chembl_id = _normalize_identifier(_text(row.get("target_chembl_id")))
                symbol = _canonical_hgnc_symbol(row.get("canonical_hgnc_symbol"))
                if not target_chembl_id or not symbol:
                    if target_chembl_id and not symbol:
                        invalid_symbols += 1
                    continue
                grouped.setdefault(target_chembl_id, []).append({
                    **dict(row),
                    "canonical_hgnc_symbol": symbol,
                })
        accepted = 0
        ambiguous = 0
        for target_chembl_id, rows in grouped.items():
            symbols = {_canonical_hgnc_symbol(row.get("canonical_hgnc_symbol")) for row in rows}
            symbols.discard("")
            component_counts = {_int(row.get("target_component_count")) for row in rows}
            component_counts.discard(None)
            if len(symbols) != 1 or component_counts != {1}:
                ambiguous += 1
                continue
            row = rows[0]
            self.mapping[target_chembl_id] = {
                "target_chembl_id": _text(row.get("target_chembl_id")),
                "target_name": _text(row.get("target_name")),
                "target_symbol": next(iter(symbols)),
                "target_symbol_lookup_key": _gene_symbol_lookup_key(next(iter(symbols))),
                "hgnc_id": _text(row.get("hgnc_id")),
                "uniprot_id": _text(row.get("component_uniprot_id")),
                "chembl_target_type": _text(row.get("chembl_target_type")),
                "component_id": _text(row.get("component_id")),
                "mapping_method": _text(row.get("mapping_method")),
                "mapping_status": _text(row.get("mapping_status")),
                "mapping_reason": _text(row.get("mapping_reason")),
                "source": self.source,
                "release": self.release,
            }
            accepted += 1
        self.status = "loaded"
        self.load_warnings.extend([
            f"chembl_target_mapping_ambiguous_or_multicomponent:{ambiguous}",
            f"chembl_target_mapping_invalid_hgnc_symbol:{invalid_symbols}",
        ])
        if accepted == 0:
            self.load_warnings.append("chembl_target_mapping_no_single_component_targets")

    def enrich(
        self,
        targets: Sequence[PlatoTargetEvidence],
    ) -> Tuple[List[PlatoTargetEvidence], Dict[str, Any]]:
        if self.status != "loaded":
            return list(targets), {
                "status": self.status,
                "path": str(self.path) if self.path else "",
                "sha256": self.sha256,
                "mapped_targets": 0,
                "warnings": list(self.load_warnings),
            }
        enriched: List[PlatoTargetEvidence] = []
        mapped = 0
        unmapped = 0
        for target in targets:
            clone = PlatoTargetEvidence(**asdict(target))
            key = _normalize_identifier(clone.target_chembl_id)
            mapping = self.mapping.get(key)
            if mapping:
                if not clone.target_symbol:
                    clone.target_symbol = mapping["target_symbol"]
                if not clone.uniprot_id:
                    clone.uniprot_id = mapping.get("uniprot_id", "")
                clone.raw = {
                    **dict(clone.raw or {}),
                    "chembl_target_mapping": mapping,
                }
                if "chembl_r37_target_mapping_applied" not in clone.parse_warnings:
                    clone.parse_warnings.append("chembl_r37_target_mapping_applied")
                mapped += 1
            else:
                clone.raw = {
                    **dict(clone.raw or {}),
                    "chembl_target_mapping": {
                        "status": "not_mapped_or_multicomponent",
                        "release": self.release,
                    },
                }
                unmapped += 1
            enriched.append(clone)
        return enriched, {
            "status": "loaded",
            "path": str(self.path) if self.path else "",
            "source": self.source,
            "release": self.release,
            "sha256": self.sha256,
            "single_component_mapping_count": len(self.mapping),
            "mapped_targets": mapped,
            "unmapped_or_multicomponent_targets": unmapped,
            "warnings": list(self.load_warnings),
            "semantics": "Identifier translation only; target candidacy remains from PLATO.",
        }


def _stable_json(value: Any) -> str:
    return json.dumps(
        _json_safe_value(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _plato_duplicate_keys(target: PlatoTargetEvidence) -> List[str]:
    """Return identity keys used only to detect repeated provider assertions."""
    keys: List[str] = []
    chembl = _normalize_identifier(target.target_chembl_id)
    uniprot = _normalize_identifier(target.uniprot_id)
    # Keep punctuation in the duplicate key so a malformed ``HLAA`` assertion
    # cannot silently overwrite the approved ``HLA-A`` spelling.
    symbol = _canonical_hgnc_symbol(target.target_symbol)
    if chembl:
        keys.append(f"chembl:{chembl}")
    if uniprot:
        keys.append(f"uniprot:{uniprot}")
    if symbol:
        keys.append(f"symbol:{symbol}")
    if not keys:
        name = _normalize_name(target.target_name)
        organism = _normalize_name(target.organism)
        taxid = _taxid(target.taxid)
        if name:
            keys.append(f"fallback:{name}|{organism}|{taxid}")
    return keys


def _deterministic_target_key(target: PlatoTargetEvidence) -> Tuple[Any, ...]:
    score = target.score if target.score is not None else float("-inf")
    rank_is_explicit = "rank_inferred_from_payload_order" not in target.parse_warnings
    rank = target.rank if rank_is_explicit and target.rank is not None else 2**31 - 1
    identity = (
        _normalize_identifier(target.target_chembl_id),
        _normalize_identifier(target.uniprot_id),
        _canonical_hgnc_symbol(target.target_symbol),
        _normalize_name(target.target_name),
        _normalize_name(target.organism),
        _taxid(target.taxid),
    )
    return (-score, rank, identity, _stable_json(target.raw))


def _consistent_duplicate_value(
    group: Sequence[PlatoTargetEvidence],
    attribute: str,
    normalize,
) -> Tuple[str, bool]:
    observed: Dict[str, List[str]] = {}
    for target in group:
        raw_value = _text(getattr(target, attribute))
        if not raw_value:
            continue
        normalized = normalize(raw_value)
        if normalized:
            observed.setdefault(normalized, []).append(raw_value)
    if len(observed) != 1:
        return "", len(observed) > 1
    # Equivalent spellings are resolved without relying on payload order.
    values = next(iter(observed.values()))
    return sorted(set(values), key=lambda item: (item.casefold(), item))[0], False


def _merge_duplicate_plato_group(
    group: Sequence[PlatoTargetEvidence],
) -> PlatoTargetEvidence:
    """Merge repeated assertions deterministically and conservatively.

    Conflicting identity/species/type assertions are cleared.  Downstream
    human/protein filters will therefore fail closed instead of accepting the
    first JSON record encountered.
    """
    if len(group) == 1:
        return group[0]
    warnings = sorted({warning for target in group for warning in target.parse_warnings})
    warnings.extend((
        "duplicate_plato_records_merged_deterministically",
        f"duplicate_plato_record_count:{len(group)}",
    ))

    fields = {
        "target_chembl_id": lambda value: _normalize_identifier(value),
        "uniprot_id": lambda value: _normalize_identifier(value),
        "target_symbol": lambda value: _canonical_hgnc_symbol(value),
        "organism": lambda value: _normalize_name(value),
        "taxid": lambda value: _taxid(value),
        "target_type": lambda value: _normalize_name(value),
        "action": lambda value: _normalize_name(value),
    }
    merged_values: Dict[str, str] = {}
    conflicting_fields: set[str] = set()
    for attribute, normalizer in fields.items():
        merged_values[attribute], conflict = _consistent_duplicate_value(
            group, attribute, normalizer,
        )
        if conflict:
            conflicting_fields.add(attribute)
            warnings.append(f"duplicate_plato_conflicting_{attribute}_cleared")
    if conflicting_fields & {"target_chembl_id", "uniprot_id", "target_symbol"}:
        # Do not allow a target-name fallback to turn contradictory strong
        # identifiers into mechanism evidence.
        merged_values["target_type"] = ""
        warnings.append("duplicate_plato_identity_conflict_failed_closed")
    if conflicting_fields & {"organism", "taxid"}:
        merged_values["organism"] = ""
        merged_values["taxid"] = ""
        warnings.append("duplicate_plato_species_conflict_failed_closed")
    if merged_values["target_symbol"]:
        merged_values["target_symbol"] = _canonical_hgnc_symbol(
            merged_values["target_symbol"]
        )

    names = sorted(
        {_text(target.target_name) for target in group if _text(target.target_name)},
        key=lambda item: (_normalize_name(item), item.casefold(), item),
    )
    target_name = names[0] if names else ""
    if len({_normalize_name(name) for name in names}) > 1:
        warnings.append("duplicate_plato_target_name_disagreement_deterministic_choice")

    scores = [target.score for target in group]
    score = min(scores) if all(value is not None for value in scores) else None
    if len(set(scores)) > 1:
        warnings.append("duplicate_plato_score_disagreement_conservative_minimum")
    reliable_values = [target.reliable for target in group]
    reliable = reliable_values[0] if len(set(reliable_values)) == 1 else None
    if len(set(reliable_values)) > 1:
        warnings.append("duplicate_plato_reliability_disagreement_cleared")
    reliable_raw_values = sorted({
        _text(target.reliable_raw) for target in group if _text(target.reliable_raw)
    }, key=lambda item: (item.casefold(), item))
    reliable_raw = reliable_raw_values[0] if len(reliable_raw_values) == 1 else ""

    applicability_values = [target.applicability_score for target in group]
    applicability = (
        min(applicability_values)
        if all(value is not None for value in applicability_values)
        else None
    )
    if len(set(applicability_values)) > 1:
        warnings.append("duplicate_plato_applicability_disagreement_conservative_minimum")

    explicit_ranks = [
        target.rank for target in group
        if target.rank is not None
        and "rank_inferred_from_payload_order" not in target.parse_warnings
    ]
    all_ranks_explicit = len(explicit_ranks) == len(group)
    rank = max(explicit_ranks) if all_ranks_explicit else None
    if not all_ranks_explicit:
        warnings.append("duplicate_plato_rank_discarded_as_order_dependent")
    elif len(set(explicit_ranks)) > 1:
        warnings.append("duplicate_plato_rank_disagreement_conservative_worst")

    ligand_map: Dict[str, Dict[str, Any]] = {}
    for target in group:
        for ligand in target.supporting_ligands:
            safe_ligand = _json_safe_value(ligand)
            if isinstance(safe_ligand, Mapping):
                ligand_map.setdefault(_stable_json(safe_ligand), dict(safe_ligand))

    raw_records = sorted(
        (_json_safe_value(target.raw) for target in group),
        key=_stable_json,
    )
    return PlatoTargetEvidence(
        target_name=(
            target_name
            or merged_values["target_symbol"]
            or merged_values["target_chembl_id"]
            or merged_values["uniprot_id"]
        ),
        target_symbol=merged_values["target_symbol"],
        target_chembl_id=merged_values["target_chembl_id"],
        uniprot_id=merged_values["uniprot_id"],
        organism=merged_values["organism"],
        target_type=merged_values["target_type"],
        taxid=merged_values["taxid"],
        score=score,
        reliable=reliable,
        reliable_raw=reliable_raw,
        rank=rank,
        applicability_score=applicability,
        action=merged_values["action"],
        supporting_ligands=[ligand_map[key] for key in sorted(ligand_map)],
        parse_warnings=list(dict.fromkeys(warnings)),
        raw={
            "duplicate_merge_policy": "deterministic_conservative_v1",
            "duplicate_records": raw_records,
        },
    )


def _deduplicate_plato_targets(
    targets: Sequence[PlatoTargetEvidence],
) -> List[PlatoTargetEvidence]:
    """Coalesce records sharing a strong identifier using order-free groups."""
    if not targets:
        return []
    parents = list(range(len(targets)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    key_owner: Dict[str, int] = {}
    for index, target in enumerate(targets):
        for key in _plato_duplicate_keys(target):
            owner = key_owner.setdefault(key, index)
            union(index, owner)

    grouped: Dict[int, List[PlatoTargetEvidence]] = {}
    for index, target in enumerate(targets):
        grouped.setdefault(find(index), []).append(target)
    merged = [_merge_duplicate_plato_group(group) for group in grouped.values()]
    return sorted(merged, key=_deterministic_target_key)


def _candidate_records(payload: Any) -> Iterable[Mapping[str, Any]]:
    """Yield target-like records from PLATO or normalized JSON layouts."""
    if isinstance(payload, list):
        for item in payload:
            yield from _candidate_records(item)
        return
    if not isinstance(payload, dict):
        return

    target_keys = {
        "target", "target_name", "name", "pref_name", "target_chembl_id",
        "chembl_target_id", "target_id", "uniprot", "uniprot_id",
        "chemblid", "nome",
    }
    lowered = {str(key).strip().lower() for key in payload}
    if lowered & target_keys and not lowered <= {"name", "value"}:
        yield payload

    for key, value in payload.items():
        if str(key).strip().lower() in {
            "targets", "predictions", "results", "target_fishing", "targetfishing",
            "data", "hits", "response",
        }:
            yield from _candidate_records(value)


def parse_plato_json(payload: Any) -> List[PlatoTargetEvidence]:
    """Parse PLATO JSON while tolerating release-specific field aliases."""
    parsed: List[PlatoTargetEvidence] = []
    for record_index, record in enumerate(_candidate_records(payload), start=1):
        raw_record = dict(record)
        nested_target = _first(record, ("target",), None)
        if isinstance(nested_target, Mapping):
            record = {**dict(nested_target), **dict(record)}
        parse_warnings: List[str] = []
        target_name = _label(_first(record, ("target_name", "target", "pref_name", "name", "nome")))
        chembl_id = _text(_first(record, ("target_chembl_id", "chembl_target_id", "chembl_id", "chemblid")))
        if "nome" in {str(key).strip().lower() for key in record}:
            parts = target_name.rsplit(":", 1)
            if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                target_name = _text(parts[0])
                record = {**dict(record), "organism": _text(parts[1])}
                parse_warnings.append("organism_parsed_from_plato_r37_nome")
        raw_target_id = _text(_first(record, ("target_id",)))
        if raw_target_id.upper().startswith("CHEMBL") and not chembl_id:
            chembl_id = raw_target_id
            parse_warnings.append("target_id_interpreted_as_chembl_id")
        symbol = _text(_first(record, ("target_symbol", "gene_symbol", "gene", "symbol")))
        uniprot = _text(_first(record, ("uniprot_id", "uniprot", "accession")))
        if raw_target_id and not chembl_id and not uniprot and _looks_like_uniprot_accession(raw_target_id):
            uniprot = raw_target_id
            parse_warnings.append("target_id_interpreted_as_uniprot")
        elif raw_target_id and not chembl_id and not symbol and _looks_like_gene_symbol(raw_target_id):
            symbol = raw_target_id
            parse_warnings.append("target_id_interpreted_as_gene_symbol")
        elif raw_target_id and not any((target_name, chembl_id, symbol, uniprot)):
            target_name = raw_target_id
            parse_warnings.append("untyped_target_id_preserved_as_name")
        organism = _label(_first(record, ("organism", "species", "target_organism"), ""))
        taxid = _taxid(_first(record, ("taxid", "tax_id", "taxonomy_id", "organism_id"), ""))
        target_type = _text(_first(record, ("target_type", "type")))
        if not target_type and _text(_first(record, ("chemblid",))) and _text(_first(record, ("nome",))):
            target_type = "single protein"
            parse_warnings.append("target_type_assigned_from_plato_r37_target_fishing_schema")
        if not organism and not taxid:
            parse_warnings.append("missing_organism")
        if not target_type:
            parse_warnings.append("missing_target_type")
        else:
            normalized_target_type = _normalize_name(target_type)
            if normalized_target_type in NON_PROTEIN_TARGET_TYPES or re.search(
                r"\b(?:non|not)(?:\s+a)?\s+protein\b|\bnonprotein\b", normalized_target_type,
            ):
                parse_warnings.append("explicit_non_protein_target_type")
            elif normalized_target_type not in PROTEIN_TARGET_TYPES:
                parse_warnings.append("unrecognized_protein_target_type")
        raw_score = _first(record, ("score", "plato_score", "target_score", "similarity_score"), None)
        score, score_warnings = _validated_plato_score(raw_score)
        parse_warnings.extend(score_warnings)
        reliable_raw_value = _first(record, ("reliable", "reliability", "is_reliable", "validato"), None)
        reliable = _bool(reliable_raw_value)
        if reliable_raw_value in (None, ""):
            parse_warnings.extend(("reliability_missing", "reliability_unknown"))
        elif reliable is None:
            parse_warnings.extend(("reliability_unrecognized", "reliability_unknown"))
        elif reliable is False:
            parse_warnings.append("reliability_false")
        if organism and taxid:
            taxid_human = _taxid(taxid) == "9606"
            if taxid_human != _is_human_label(organism):
                parse_warnings.append("species_conflict")
        elif organism and not _is_human_label(organism):
            parse_warnings.append("non_human_or_uncontrolled_organism_label")
        raw_rank = _first(record, ("rank", "ranking", "position"), None)
        rank = _int(raw_rank) if raw_rank not in (None, "") else record_index
        if raw_rank in (None, ""):
            parse_warnings.append("rank_inferred_from_payload_order")
        applicability = _float(_first(record, ("max_similarity", "similarity", "tanimoto", "applicability_score"), None))
        action = _text(_first(record, ("action", "action_type", "mechanism_of_action", "moa")))
        ligands = _first(record, ("supporting_ligands", "similar_ligands", "ligands", "compounds", "simili"), [])
        if not isinstance(ligands, list):
            ligands = []

        if not any((target_name, chembl_id, symbol, uniprot)):
            continue
        parsed.append(PlatoTargetEvidence(
            target_name=target_name or symbol or chembl_id or uniprot,
            target_symbol=symbol,
            target_chembl_id=chembl_id,
            uniprot_id=uniprot,
            organism=organism,
            target_type=target_type,
            taxid=taxid,
            score=score,
            reliable=reliable,
            reliable_raw=_text(reliable_raw_value),
            rank=rank,
            applicability_score=applicability,
            action=action,
            supporting_ligands=[dict(item) for item in ligands if isinstance(item, dict)],
            parse_warnings=parse_warnings,
            raw=raw_record,
        ))
    return _deduplicate_plato_targets(parsed)


def _payload_header_value(payload: Any, names: Sequence[str]) -> str:
    """Read contract metadata only from the payload header, never target rows."""
    if not isinstance(payload, Mapping):
        return ""
    value = _label(_first(payload, names, ""))
    if value:
        return value
    metadata = _first(payload, ("metadata", "meta", "header"), None)
    return _label(_first(metadata, names, "")) if isinstance(metadata, Mapping) else ""


def _json_schema_fingerprint(payload: Any) -> Tuple[str, List[str]]:
    """Fingerprint JSON key paths and value kinds without using scalar values."""
    signatures: set[str] = set()

    def visit(value: Any, path: str, depth: int = 0) -> None:
        if depth > 16:
            signatures.add(f"{path}:depth_limit")
            return
        if isinstance(value, Mapping):
            signatures.add(f"{path}:object")
            for key in sorted(value, key=lambda item: str(item).casefold()):
                normalized_key = str(key).strip().lower()
                visit(value[key], f"{path}.{normalized_key}", depth + 1)
            return
        if isinstance(value, list):
            signatures.add(f"{path}:array")
            for item in value:
                visit(item, f"{path}[]", depth + 1)
            return
        if value is None:
            kind = "null"
        elif isinstance(value, bool):
            kind = "boolean"
        elif isinstance(value, (int, float)):
            kind = "number"
        elif isinstance(value, str):
            kind = "string"
        else:
            kind = type(value).__name__.lower()
        signatures.add(f"{path}:{kind}")

    visit(payload, "$")
    ordered = sorted(signatures)
    return sha256("\n".join(ordered).encode("utf-8")).hexdigest(), ordered


def inspect_plato_payload(payload: Any, parsed: Optional[Sequence[PlatoTargetEvidence]] = None) -> Dict[str, Any]:
    """Describe the parsed contract without claiming an unverified r37 schema."""
    records = list(_candidate_records(payload))
    targets = list(parsed) if parsed is not None else parse_plato_json(payload)
    fingerprint, schema_signature = _json_schema_fingerprint(payload)
    fingerprint_allowlisted = fingerprint in VERIFIED_PLATO_SCHEMA_FINGERPRINTS
    verified_plato_r37_target_fishing = bool(
        fingerprint_allowlisted and _looks_like_plato_r37_target_fishing_payload(payload)
    )
    reported_provider = _payload_header_value(
        payload, ("provider", "provider_name", "platform", "service"),
    )
    reported_release = _payload_header_value(payload, ("release", "version", "plato_release"))
    method = _payload_header_value(payload, ("method", "prediction_method", "mode"))
    release_source = "payload" if reported_release else ""
    if verified_plato_r37_target_fishing:
        reported_provider = reported_provider or "PLATO"
        reported_release = reported_release or "r37"
        method = method or "target_fishing"
        release_source = release_source or "verified_schema_fixture"
    normalized_provider = _normalize_name(reported_provider)
    normalized_method = _normalize_name(method)
    provider_validation_status = (
        "valid_plato" if normalized_provider in PLATO_PROVIDER_LABELS
        else "missing" if not normalized_provider else "invalid_provider"
    )
    method_validation_status = (
        "valid_target_fishing" if normalized_method in PLATO_TARGET_FISHING_METHODS
        else "missing" if not normalized_method else "invalid_method"
    )
    field_names = sorted({str(key).strip().lower() for record in records for key in record})
    schema_fixture_verified = bool(
        fingerprint_allowlisted
        and provider_validation_status == "valid_plato"
        and method_validation_status == "valid_target_fishing"
    )
    adapter = "plato_verified_fixture_adapter" if schema_fixture_verified else "generic_alias_adapter"
    warnings: List[str] = []
    if provider_validation_status == "missing":
        warnings.append("provider_not_reported_by_payload")
    elif provider_validation_status != "valid_plato":
        warnings.append("payload_provider_is_not_plato")
    if not reported_release:
        warnings.append("provider_release_not_reported_by_payload")
    if method_validation_status == "missing":
        warnings.append("prediction_method_not_reported_by_payload")
    elif method_validation_status != "valid_target_fishing":
        warnings.append("prediction_method_is_not_target_fishing")
    if not fingerprint_allowlisted:
        warnings.append("schema_fingerprint_not_verified_against_fixture")
    if any(target.human_status == "unknown_species" for target in targets):
        warnings.append("records_with_unknown_species")
    if any(not target.target_type for target in targets):
        warnings.append("records_with_unknown_target_type")
    if any(target.applicability_score is None for target in targets):
        warnings.append("applicability_domain_not_reported")
    record_warning_counts: Dict[str, int] = {}
    for target in targets:
        for warning in target.parse_warnings:
            record_warning_counts[warning] = record_warning_counts.get(warning, 0) + 1
    warnings.extend(f"target_record_warning:{warning}" for warning in sorted(record_warning_counts))
    return {
        "parser_version": PARSER_VERSION,
        "adapter": adapter,
        "provider_name": reported_provider or "unknown",
        "provider_validation_status": provider_validation_status,
        "provider_method": method or "unknown",
        "method_validation_status": method_validation_status,
        "schema_fingerprint": fingerprint,
        "schema_fingerprint_algorithm": "sha256_json_key_paths_and_value_kinds_v1",
        "schema_signature": schema_signature,
        "schema_fingerprint_allowlisted": fingerprint_allowlisted,
        "schema_fixture_verified": schema_fixture_verified,
        "reported_release": reported_release,
        "reported_release_source": release_source,
        "method": method,
        "raw_candidate_count": len(records),
        "parsed_count": len(targets),
        "field_names": field_names,
        "field_coverage": {
            "organism": sum(bool(target.organism or target.taxid) for target in targets) / max(1, len(targets)),
            "target_type": sum(bool(target.target_type) for target in targets) / max(1, len(targets)),
            "score": sum(target.score is not None for target in targets) / max(1, len(targets)),
            "reliability": sum(target.reliable is not None for target in targets) / max(1, len(targets)),
            "applicability": sum(target.applicability_score is not None for target in targets) / max(1, len(targets)),
        },
        "human_status_counts": {
            status: sum(target.human_status == status for target in targets)
            for status in sorted({target.human_status for target in targets})
        },
        "record_warning_counts": record_warning_counts,
        "validation": {
            "provider": provider_validation_status,
            "method": method_validation_status,
            "schema": "verified_fixture" if schema_fixture_verified else "unverified_generic_schema",
            "platform_authenticity": "not_authenticated",
        },
        "warnings": warnings,
    }


def _atomic_write_new(path: Path, data: bytes) -> bool:
    """Atomically create an immutable file; never replace different content."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() == data:
            return False
        raise FileExistsError(f"Refusing to overwrite immutable cache artifact: {path}")

    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            # A hard-link publish is atomic and fails if another writer won.
            os.link(temporary, path)
        except FileExistsError:
            if path.read_bytes() != data:
                raise FileExistsError(f"Concurrent immutable cache conflict: {path}")
            return False
        except OSError:
            # Preserve no-overwrite semantics on filesystems without links.
            if path.exists():
                if path.read_bytes() != data:
                    raise FileExistsError(f"Concurrent immutable cache conflict: {path}")
                return False
            try:
                with path.open("xb") as handle:
                    handle.write(data)
                    handle.flush()
                    os.fsync(handle.fileno())
            except FileExistsError:
                if path.read_bytes() != data:
                    raise FileExistsError(f"Concurrent immutable cache conflict: {path}")
                return False
        return True
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


class PlatoPredictionCache:
    """Filesystem cache for raw PLATO JSON and normalized metadata."""

    def __init__(self, root: str | Path):
        self.root = Path(root)

    def path_for(self, smiles: str) -> Path:
        return self.root / f"{smiles_cache_key(smiles)}.json"

    def _candidate_paths(self, smiles: str) -> List[Path]:
        canonical = self.path_for(smiles)
        legacy = self.root / f"{_legacy_smiles_cache_key(smiles)}.json"
        return [canonical] if canonical == legacy else [canonical, legacy]

    def load_ignoring_stereochemistry(self, smiles: str) -> Tuple[List[PlatoTargetEvidence], Dict[str, Any]]:
        """Reuse a consistent connectivity group, retaining exact query provenance.

        Original cache keys and raw records remain stereospecific and immutable.
        Conflicting predictions in a group are returned as an explicit conflict.
        """
        requested = standardize_smiles(smiles)
        group = stereo_agnostic_smiles(smiles)
        exact_targets, exact = self.load(smiles)
        matching = {"policy": "ignore_stereochemistry", "requested_structure": requested,
                    "canonical_group_smiles": group, "isotopes_charge_and_components_preserved": True}
        candidates: Dict[str, str] = {}
        if str(exact.get("status", "")).startswith("cache_hit"):
            # Never hide a corrupted exact entry behind a different group member.
            exact_query = exact.get("query", {})
            if (not structure_identity_matches(exact_query, requested)
                    or exact.get("parser_warnings") and any(
                        "mismatch" in str(value) or "integrity" in str(value)
                        for value in exact["parser_warnings"])):
                exact["structure_matching"] = {**matching, "status": "invalid_exact_entry"}
                return exact_targets, exact
            candidates[exact_query["canonical_isomeric_smiles"]] = exact["cache_path"]
        for path in sorted(self.root.glob("*.json")):
            try:
                envelope = _loads_strict_json(path.read_text(encoding="utf-8"))
            except (OSError, ValueError, UnicodeError):
                continue  # Unrelated exports are not cache entries.
            if not isinstance(envelope, Mapping) or envelope.get("alias_to"):
                continue
            query = envelope.get("query", {})
            if not isinstance(query, Mapping):
                continue
            source_smiles = _text(query.get("canonical_isomeric_smiles"))
            if not source_smiles or not envelope.get("cache_schema_version"):
                continue
            try:
                matches = stereo_agnostic_smiles(source_smiles) == group
            except (ValueError, ImportError):
                matches = False
            if matches and path.resolve() in {candidate.resolve() for candidate in self._candidate_paths(source_smiles)}:
                candidates[source_smiles] = str(path)
        if not candidates:
            exact["structure_matching"] = {**matching, "status": "no_group_cache"}
            return exact_targets, exact
        members, loaded = [], []
        for source_smiles in sorted(candidates, key=lambda value: (value != group, value)):
            parsed, metadata = self.load(source_smiles)
            query = metadata.get("query", {})
            # Rank records provider order; compare substantive target evidence.
            snapshots = [{key: value for key, value in target.to_dict().items() if key != "rank"}
                         for target in parsed]
            fingerprint = _payload_sha256(sorted(snapshots, key=lambda value: json.dumps(value, sort_keys=True)))
            member = {"cache_path": metadata.get("cache_path"), "query_smiles": source_smiles,
                      "inchikey": query.get("inchikey"), "raw_sha256": metadata.get("source", {}).get("sha256"),
                      "target_snapshot_sha256": fingerprint}
            members.append(member)
            if (not structure_identity_matches(query, standardize_smiles(source_smiles))
                    or metadata.get("raw_blob_integrity_status") != "verified_immutable_raw_blob"
                    or any("mismatch" in str(value) or "integrity" in str(value)
                           for value in metadata.get("parser_warnings", []))):
                metadata["structure_matching"] = {**matching, "status": "invalid_group_member", "members": members}
                metadata["status"] = "cache_error_structure_group_integrity"
                return [], metadata
            loaded.append((parsed, metadata))
        if len({member["target_snapshot_sha256"] for member in members}) != 1:
            return [], {"status": "cache_error_stereoisomer_prediction_conflict",
                        "error": "Same connectivity has different PLATO target results; no arbitrary merge performed",
                        "structure_matching": {**matching, "status": "conflicting_predictions", "members": members}}
        parsed, metadata = loaded[0]
        metadata["structure_matching"] = {
            **matching, "status": "matched", "members": members, "n_cached_structures": len(members),
            "source_query_structure": dict(metadata["query"]),
            "cross_stereochemistry_reuse": not structure_identity_matches(metadata["query"], requested),
            "consistent_cached_target_results": True,
            "equivalent_biological_activity_claimed": False,
        }
        return parsed, metadata

    def load(self, smiles: str) -> Tuple[List[PlatoTargetEvidence], Dict[str, Any]]:
        paths = self._candidate_paths(smiles)
        path = next((candidate for candidate in paths if candidate.exists()), paths[0])
        if not path.exists():
            return [], {
                "status": "cache_miss",
                "cache_path": str(path),
                "candidate_paths": [str(candidate) for candidate in paths],
                "query": standardize_smiles(smiles),
            }
        envelope = _loads_strict_json(path.read_text(encoding="utf-8"))
        if not isinstance(envelope, Mapping):
            raise ValueError(f"PLATO cache envelope must be an object: {path}")
        alias_path: Optional[Path] = None
        if envelope.get("alias_to"):
            alias_path = path
            alias_value = envelope.get("alias_to")
            if not isinstance(alias_value, str) or not alias_value.strip():
                raise ValueError(f"PLATO cache alias target must be a filename string: {path}")
            target_name = Path(alias_value.strip()).name
            if target_name != alias_value.strip() or not target_name.lower().endswith(".json"):
                raise ValueError(f"PLATO cache alias target must be a local JSON filename: {path}")
            path = self.root / target_name
            root_resolved = self.root.resolve()
            target_resolved = path.resolve()
            if target_resolved != root_resolved and root_resolved not in target_resolved.parents:
                raise ValueError(f"PLATO cache alias target resolves outside cache root: {path}")
            if not path.exists():
                raise ValueError(f"PLATO cache alias target is missing: {path}")
            envelope = _loads_strict_json(path.read_text(encoding="utf-8"))
            if not isinstance(envelope, Mapping):
                raise ValueError(f"PLATO cache alias target must be an object: {path}")
            if envelope.get("alias_to"):
                raise ValueError(f"Nested PLATO cache aliases are not allowed: {path}")
        cache_schema_version = _text(envelope.get("cache_schema_version", "1"))
        if cache_schema_version == CACHE_SCHEMA_VERSION and "raw_result" not in envelope:
            raise ValueError(f"Current PLATO cache envelope is missing raw_result: {path}")
        raw_result = envelope.get("raw_result", envelope)
        normalized_targets = envelope.get("normalized_targets")
        targets = parse_plato_json(raw_result)
        runtime_report = inspect_plato_payload(raw_result, targets)
        warnings = list(runtime_report.get("warnings", []))
        metadata = _cache_mapping(envelope.get("metadata"), "metadata", warnings)
        provider = _cache_mapping(envelope.get("provider"), "provider", warnings)
        query = _cache_mapping(envelope.get("query"), "query", warnings)
        source = _cache_mapping(envelope.get("source"), "source", warnings)
        parser = _cache_mapping(envelope.get("parser"), "parser", warnings)
        stored_warnings = _cache_string_list(parser.get("warnings"), "parser_warnings", warnings)
        for field in (
            "adapter", "provider_name", "provider_validation_status", "provider_method",
            "method_validation_status", "schema_fingerprint", "schema_fixture_verified",
            "reported_release", "method", "raw_candidate_count", "parsed_count",
            "field_names", "warnings",
        ):
            if field in parser and parser.get(field) != runtime_report.get(field):
                warnings.append(f"cached_parser_report_{field}_mismatch")
        expected_payload_hash = _text(source.get("payload_sha256")).lower()
        if expected_payload_hash and not _is_sha256(expected_payload_hash):
            warnings.append("cached_payload_hash_invalid_format")
            warnings.append("cached_payload_hash_mismatch")
        elif expected_payload_hash and expected_payload_hash != _payload_sha256(raw_result):
            warnings.append("cached_payload_hash_mismatch")
        if cache_schema_version != CACHE_SCHEMA_VERSION:
            warnings.append("cache_schema_not_current")
        raw_blob_integrity_status = "not_available_legacy_cache"
        raw_blob_relative_path = _text(source.get("raw_blob_relative_path"))
        if cache_schema_version == CACHE_SCHEMA_VERSION:
            raw_blob_integrity_status = "raw_blob_integrity_failure"
            if not raw_blob_relative_path:
                warnings.append("raw_blob_path_missing")
            else:
                try:
                    root_resolved = self.root.resolve()
                    raw_blob_path = (self.root / raw_blob_relative_path).resolve()
                    if raw_blob_path != root_resolved and root_resolved not in raw_blob_path.parents:
                        warnings.append("raw_blob_path_outside_cache_root")
                    elif not raw_blob_path.is_file():
                        warnings.append("raw_blob_missing")
                    else:
                        raw_bytes = raw_blob_path.read_bytes()
                        raw_hash = sha256(raw_bytes).hexdigest()
                        expected_raw_hash = _text(source.get("sha256"))
                        if not expected_raw_hash:
                            warnings.append("raw_blob_hash_missing")
                        elif not _is_sha256(expected_raw_hash):
                            warnings.append("raw_blob_hash_invalid_format")
                        elif raw_hash != expected_raw_hash:
                            warnings.append("raw_blob_hash_mismatch")
                        if _int(source.get("size_bytes")) != len(raw_bytes):
                            warnings.append("raw_blob_size_mismatch")
                        try:
                            raw_blob_payload = _loads_strict_json(raw_bytes.decode("utf-8-sig"))
                        except (UnicodeDecodeError, ValueError, TypeError, json.JSONDecodeError):
                            warnings.append("raw_blob_invalid_json")
                        else:
                            if _payload_sha256(raw_blob_payload) != _payload_sha256(raw_result):
                                warnings.append("raw_blob_payload_mismatch")
                        raw_blob_failures = {
                            "raw_blob_path_outside_cache_root", "raw_blob_missing", "raw_blob_hash_missing",
                            "raw_blob_hash_invalid_format", "raw_blob_hash_mismatch", "raw_blob_size_mismatch", "raw_blob_invalid_json",
                            "raw_blob_payload_mismatch",
                        }
                        if not any(warning in raw_blob_failures for warning in warnings):
                            raw_blob_integrity_status = "verified_immutable_raw_blob"
                except OSError:
                    warnings.append("raw_blob_unreadable")
        expected_normalized_hash = _text(source.get("normalized_targets_sha256")).lower()
        if cache_schema_version == CACHE_SCHEMA_VERSION:
            if not isinstance(normalized_targets, list) or not expected_normalized_hash:
                warnings.append("normalized_targets_hash_missing")
            elif not _is_sha256(expected_normalized_hash):
                warnings.append("normalized_targets_hash_invalid_format")
                warnings.append("normalized_targets_hash_mismatch")
            elif expected_normalized_hash != _payload_sha256(normalized_targets):
                warnings.append("normalized_targets_hash_mismatch")
            elif expected_normalized_hash != _payload_sha256([target.to_dict() for target in targets]):
                warnings.append("parser_output_changed_since_import")
            expected_audit_hash = _text(source.get("audit_metadata_sha256")).lower()
            audit_snapshot = {"provider": provider, "query": query, "parser": parser}
            if not expected_audit_hash:
                warnings.append("audit_metadata_hash_missing")
            elif not _is_sha256(expected_audit_hash):
                warnings.append("audit_metadata_hash_invalid_format")
                warnings.append("audit_metadata_hash_mismatch")
            elif expected_audit_hash != _payload_sha256(audit_snapshot):
                warnings.append("audit_metadata_hash_mismatch")
        if parser.get("parser_version") and parser.get("parser_version") != PARSER_VERSION:
            warnings.append("cache_parser_version_differs_from_runtime")
        requested = standardize_smiles(smiles)
        cached_canonical = _text(query.get("canonical_isomeric_smiles"))
        if cached_canonical and cached_canonical != requested["canonical_isomeric_smiles"]:
            if (
                requested.get("valid") is not True
                and _text(query.get("submitted_smiles")) == _text(smiles)
            ):
                warnings.append("canonical_identity_not_reverified_without_rdkit")
            else:
                warnings.append("cache_query_structure_mismatch")
        binding_status = _text(query.get("query_binding_status")) or "legacy_cache_missing_binding_status"
        if binding_status not in VERIFIED_QUERY_BINDING_STATUSES:
            warnings.append("plato_query_binding_user_asserted_unverified")
        try:
            runtime_binding = _query_binding_report(raw_result, smiles)
            if (
                binding_status == "payload_query_structure_verified"
                and runtime_binding.get("query_binding_status") != binding_status
            ):
                warnings.append("cached_query_binding_status_mismatch")
        except ValueError:
            warnings.append("cached_payload_query_structure_mismatch")
        runtime_release = _text(runtime_report.get("reported_release"))
        stored_release = _text(provider.get("release", metadata.get("provider_release", "unknown"))) or "unknown"
        if runtime_release and stored_release.lower() != runtime_release.lower():
            warnings.append("cached_provider_release_mismatch")
        runtime_method = _text(runtime_report.get("method"))
        stored_method = _text(provider.get("method"))
        if runtime_method and stored_method and runtime_method.lower() != stored_method.lower():
            warnings.append("cached_provider_method_mismatch")
        user_declared_raw = parser.get(
            "user_declared_real_export", parser.get("confirmed_real_export")
        )
        user_declared_real_export = _cache_bool(
            user_declared_raw, "user_declared_real_export", warnings,
        ) is True
        stored_schema_fixture_verified = _cache_bool(
            parser.get("schema_fixture_verified"), "schema_fixture_verified", warnings,
        )
        schema_fixture_verified = bool(
            runtime_report.get("schema_fixture_verified")
            and stored_schema_fixture_verified is True
            and not any(warning.startswith("cached_parser_report_") for warning in warnings)
        )
        metadata.update({
            "status": "cache_hit_alias" if alias_path else ("cache_hit_legacy_key" if path != paths[0] else "cache_hit"),
            "cache_path": str(path),
            "cache_alias_path": str(alias_path) if alias_path else "",
            "cache_schema_version": cache_schema_version,
            "parsed_targets": len(targets),
            "provider_release": runtime_release or stored_release,
            "release_source": provider.get("release_source", "legacy_metadata"),
            "provider_name": runtime_report.get("provider_name", "unknown"),
            "provider_validation_status": runtime_report.get("provider_validation_status", "missing"),
            "provider_method": runtime_report.get("provider_method", "unknown"),
            "method_validation_status": runtime_report.get("method_validation_status", "missing"),
            "schema_fingerprint": runtime_report.get("schema_fingerprint", ""),
            "adapter": runtime_report.get("adapter", "generic_alias_adapter"),
            "query": query or requested,
            "source": {**source, "raw_blob_integrity_status": raw_blob_integrity_status},
            "raw_blob_integrity_status": raw_blob_integrity_status,
            "parser": {**dict(parser), "runtime_report": runtime_report},
            "stored_parser_warnings": stored_warnings,
            "parser_warnings": list(dict.fromkeys(warnings)),
            "user_declared_real_export": user_declared_real_export,
            "confirmed_real_export": False,
            "schema_fixture_verified": schema_fixture_verified,
            "platform_authenticity_status": "not_cryptographically_authenticated",
        })
        return targets, metadata

    def import_file(
        self,
        smiles: str,
        result_file: str | Path,
        drug_name: str = "",
        provider_release: str = "unknown",
        confirmed_real_export: bool = False,
        query_binding_evidence: Optional[Mapping[str, Any]] = None,
    ) -> Path:
        source = Path(result_file)
        raw_bytes = source.read_bytes()
        raw_sha256 = sha256(raw_bytes).hexdigest()
        payload = _loads_strict_json(raw_bytes.decode("utf-8-sig"))
        parsed = parse_plato_json(payload)
        if not parsed:
            raise ValueError(f"No target records recognized in PLATO result: {source}")
        report = inspect_plato_payload(payload, parsed)
        declared_release = _text(provider_release) or "unknown"
        reported_release = _text(report.get("reported_release"))
        if reported_release and declared_release != "unknown" and reported_release.lower() != declared_release.lower():
            raise ValueError(
                f"PLATO release mismatch: payload reports {reported_release}, import requested {declared_release}"
            )
        actual_release = reported_release or declared_release
        release_source = (
            _text(report.get("reported_release_source"))
            or ("payload" if reported_release else ("cli_explicit" if declared_release != "unknown" else "unknown"))
        )
        structure = standardize_smiles(smiles)
        binding_report = _query_binding_report(payload, smiles)
        evidence_report = _single_query_binding_report(
            smiles, raw_sha256, query_binding_evidence
        )
        if (
            binding_report.get("query_binding_status") != "payload_query_structure_verified"
            and evidence_report is not None
        ):
            binding_report = evidence_report
        structure.update(binding_report)
        normalized_snapshot = [target.to_dict() for target in parsed]
        self.root.mkdir(parents=True, exist_ok=True)
        destination = self.path_for(smiles)
        imported_at = datetime.now(timezone.utc).isoformat()
        reported_provider = _text(report.get("provider_name")) or "unknown"
        reported_method = _text(report.get("provider_method")) or "unknown"
        platform_authenticity_status = (
            "user_attested_unverified" if confirmed_real_export else "not_attested"
        )
        provider_record = {
            "name": reported_provider,
            "expected_name": "PLATO",
            "validation_status": report.get("provider_validation_status", "missing"),
            "release": actual_release,
            "release_source": release_source,
            "method": reported_method,
            "expected_method": "target_fishing",
            "method_validation_status": report.get("method_validation_status", "missing"),
            "platform_authenticity_status": platform_authenticity_status,
        }
        parser_record = {
            **report,
            "user_declared_real_export": bool(confirmed_real_export),
            "schema_fixture_verified": bool(report.get("schema_fixture_verified", False)),
            "platform_authenticity_status": platform_authenticity_status,
        }
        audit_metadata_hash = _payload_sha256({
            "provider": provider_record,
            "query": structure,
            "parser": parser_record,
        })
        raw_blob_relative_path = Path("raw_blobs") / f"{raw_sha256}.json"
        raw_blob_path = self.root / raw_blob_relative_path
        _atomic_write_new(raw_blob_path, raw_bytes)
        envelope = {
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "provider": provider_record,
            "query": structure,
            "source": {
                "filename": source.name,
                "path": str(source.resolve()),
                "size_bytes": len(raw_bytes),
                "sha256": raw_sha256,
                "raw_blob_relative_path": raw_blob_relative_path.as_posix(),
                "raw_blob_storage": "content_addressed_immutable_copy",
                "payload_sha256": _payload_sha256(payload),
                "normalized_targets_sha256": _payload_sha256(normalized_snapshot),
                "audit_metadata_sha256": audit_metadata_hash,
                "imported_at_utc": imported_at,
            },
            "parser": parser_record,
            "normalized_targets": normalized_snapshot,
            "metadata": {
                "provider": reported_provider,
                "expected_provider": "PLATO",
                "provider_validation_status": report.get("provider_validation_status", "missing"),
                "provider_release": actual_release,
                "release_source": release_source,
                "method": reported_method,
                "method_validation_status": report.get("method_validation_status", "missing"),
                "schema_fixture_verified": bool(report.get("schema_fixture_verified", False)),
                "platform_authenticity_status": platform_authenticity_status,
                "drug_name": drug_name,
                "query_smiles": smiles,
                "canonical_smiles": structure["canonical_isomeric_smiles"],
                "inchikey": structure["inchikey"],
                "imported_at": imported_at,
                "source_file": str(source.resolve()),
                "parsed_targets": len(parsed),
            },
            "raw_result": payload,
        }
        envelope_bytes = json.dumps(
            envelope, ensure_ascii=False, indent=2, allow_nan=False,
        ).encode("utf-8")
        raw_alias = self.root / f"{_legacy_smiles_cache_key(smiles)}.json"
        if destination.exists():
            existing = _loads_strict_json(destination.read_text(encoding="utf-8"))
            existing_source = existing.get("source", {}) if isinstance(existing, Mapping) else {}
            if (
                isinstance(existing_source, Mapping)
                and _text(existing.get("cache_schema_version")) == CACHE_SCHEMA_VERSION
                and _text(existing_source.get("sha256")) == raw_sha256
            ):
                if raw_alias != destination and not raw_alias.exists():
                    alias_bytes = json.dumps({
                        "cache_schema_version": f"{CACHE_SCHEMA_VERSION}-alias",
                        "alias_to": destination.name,
                        "submitted_smiles": smiles,
                        "created_at_utc": imported_at,
                    }, ensure_ascii=False, indent=2).encode("utf-8")
                    _atomic_write_new(raw_alias, alias_bytes)
                return destination
            raise FileExistsError(
                "A different or legacy cache entry already exists for this structure; "
                f"refusing silent overwrite: {destination}"
            )
        if raw_alias != destination and raw_alias.exists():
            existing_alias = _loads_strict_json(raw_alias.read_text(encoding="utf-8"))
            if not isinstance(existing_alias, Mapping) or Path(
                _text(existing_alias.get("alias_to")),
            ).name != destination.name:
                raise FileExistsError(f"Conflicting immutable cache alias: {raw_alias}")

        _atomic_write_new(destination, envelope_bytes)
        if raw_alias != destination:
            alias_bytes = json.dumps({
                "cache_schema_version": f"{CACHE_SCHEMA_VERSION}-alias",
                "alias_to": destination.name,
                "submitted_smiles": smiles,
                "created_at_utc": imported_at,
            }, ensure_ascii=False, indent=2).encode("utf-8")
            _atomic_write_new(raw_alias, alias_bytes)
        return destination


def _local_record(value: Any) -> Dict[str, Any]:
    legacy_tuple = not isinstance(value, Mapping)
    if isinstance(value, Mapping):
        record = dict(value)
        target_id = _text(_first(record, ("target_id", "target_symbol", "gene_symbol", "symbol", "uniprot_id", "target_chembl_id", "chembl_id")))
        target_name = _text(_first(record, ("target_name", "name", "pref_name"), target_id))
    else:
        parts = list(value) if isinstance(value, (tuple, list)) else [value]
        target_id = _text(parts[0] if parts else "")
        target_name = _text(parts[1] if len(parts) > 1 else target_id)
        record = {"target_id": target_id, "target_name": target_name}
    record["target_id"] = target_id
    record["target_name"] = target_name or target_id
    has_namespace = any(_text(_first(record, names)) for names in (
        ("target_symbol", "gene_symbol", "symbol"),
        ("uniprot_id", "uniprot", "accession"),
        ("target_chembl_id", "chembl_target_id", "chembl_id"),
    ))
    normalized = _normalize_identifier(target_id)
    if not has_namespace and normalized.startswith("CHEMBL"):
        record.setdefault("target_chembl_id", target_id)
    elif not has_namespace and _looks_like_uniprot_accession(target_id):
        record.setdefault("uniprot_id", target_id)
    elif not has_namespace:
        record.setdefault("target_symbol", target_id)
    if legacy_tuple:
        record.setdefault("organism_status", "legacy_tuple_human_assumption")
        record.setdefault("source", "local_legacy_tuple")
    return record


def _local_alignment_keys(record: Mapping[str, Any]) -> List[str]:
    keys: List[str] = []
    values = (
        ("uniprot", _first(record, ("uniprot_id", "uniprot", "accession"))),
        ("chembl", _first(record, ("target_chembl_id", "chembl_target_id", "chembl_id"))),
        ("symbol", _first(record, ("target_symbol", "gene_symbol", "symbol"))),
    )
    for prefix, value in values:
        normalized = _normalize_identifier(_text(value))
        if normalized:
            keys.append(f"{prefix}:{normalized}")
    name = _normalize_name(_text(_first(record, ("target_name", "name", "pref_name"))))
    if name:
        keys.append(f"name:{name}")
    return keys


def _local_human_status(record: Mapping[str, Any]) -> str:
    declared = _text(record.get("human_status") or record.get("organism_status")).lower()
    taxid = _normalize_identifier(_taxid(_first(record, ("taxid", "tax_id", "taxonomy_id"))))
    organism = _normalize_name(_label(_first(record, ("organism", "species"))))
    statuses: List[str] = []
    if declared in {"human", "human_explicit", "human_mapped", "human_name_inferred"}:
        statuses.append("human")
    elif declared == "non_human":
        statuses.append("non_human")
    elif declared == "species_conflict":
        return "species_conflict"
    if taxid:
        statuses.append("human" if taxid == "9606" else "non_human")
    if organism:
        statuses.append("human" if _is_human_label(organism) else "non_human")
    if "human" in statuses and "non_human" in statuses:
        return "species_conflict"
    if statuses:
        return "human_explicit" if statuses[0] == "human" else "non_human"
    if record.get("organism_status") in {"legacy_tuple_human_assumption", "local_database_human_assumption"}:
        return "human_assumed_legacy"
    return "unknown_species"


def _strong_identifiers(record: Mapping[str, Any]) -> Dict[str, str]:
    return {
        "uniprot": _normalize_identifier(_text(_first(record, ("uniprot_id", "uniprot", "accession")))),
        "chembl": _normalize_identifier(_text(_first(record, ("target_chembl_id", "chembl_target_id", "chembl_id")))),
        "taxid": _normalize_identifier(_taxid(_first(record, ("taxid", "tax_id", "taxonomy_id")))),
    }


def _identifier_conflict(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    left_ids, right_ids = _strong_identifiers(left), _strong_identifiers(right)
    return any(left_ids[key] and right_ids[key] and left_ids[key] != right_ids[key] for key in left_ids)


def _fused_identifier_conflict(current: FusedTarget, target: PlatoTargetEvidence) -> bool:
    candidate = target.to_dict()
    if current.local is not None and _identifier_conflict(current.local, candidate):
        return True
    return current.plato is not None and _identifier_conflict(current.plato.to_dict(), candidate)


def fuse_local_and_plato_targets(
    local_targets: Sequence[Any],
    plato_targets: Sequence[PlatoTargetEvidence],
    human_only: bool = True,
    reliable_only: bool = False,
    max_plato_targets: int = 20,
) -> Tuple[List[FusedTarget], Dict[str, Any]]:
    """Align human local/PLATO evidence without collapsing ambiguous records."""
    if max_plato_targets < 0:
        raise ValueError("max_plato_targets must be non-negative")
    fused: List[FusedTarget] = []
    key_to_indices: Dict[str, set[int]] = {}
    excluded_local_non_human = 0
    excluded_local_non_protein = 0
    local_target_type_unspecified = 0
    local_target_type_unknown = 0
    local_origin_unspecified = 0
    local_identifier_conflicts = 0
    target_symbol_disagreements = 0

    for raw_local in local_targets:
        record = _local_record(raw_local)
        target_id, target_name = _text(record.get("target_id")), _text(record.get("target_name"))
        if not target_id:
            continue
        human_status = _local_human_status(record)
        record["human_status"] = human_status
        if human_only and human_status in {"non_human", "species_conflict"}:
            excluded_local_non_human += 1
            continue
        local_target_type = _normalize_name(_text(_first(record, ("target_type", "type"))))
        if local_target_type:
            if (
                local_target_type in NON_PROTEIN_TARGET_TYPES
                or any(term in local_target_type for term in ("cell line", "tissue", "organism", "subcellular"))
            ):
                excluded_local_non_protein += 1
                continue
            if not (
                local_target_type in PROTEIN_TARGET_TYPES
                or "protein" in local_target_type
                or any(term in local_target_type for term in (
                    "kinase", "enzyme", "receptor", "ion channel", "transporter", "phosphatase",
                ))
            ):
                local_target_type_unknown += 1
        else:
            local_target_type_unspecified += 1
        if human_status in {"unknown_species", "human_assumed_legacy"}:
            local_origin_unspecified += 1
        keys = _local_alignment_keys(record)
        existing_indices = {
            index for key in keys for index in key_to_indices.get(key, set())
        }
        if len(existing_indices) == 1 and not _identifier_conflict(
            fused[next(iter(existing_indices))].local or {}, record
        ):
            current = fused[next(iter(existing_indices))]
            current_symbol = _normalize_identifier(_text(_first(
                current.local or {}, ("target_symbol", "gene_symbol", "symbol")
            )))
            candidate_symbol = _normalize_identifier(_text(_first(
                record, ("target_symbol", "gene_symbol", "symbol")
            )))
            if current_symbol and candidate_symbol and current_symbol != candidate_symbol:
                target_symbol_disagreements += 1
                if current.local is not None:
                    warnings = current.local.setdefault("alignment_warnings", [])
                    warning = "local_symbol_disagreement_on_shared_identifier"
                    if warning not in warnings:
                        warnings.append(warning)
                    conflicts = current.local.setdefault("symbol_conflicts", [])
                    for value in (current_symbol, candidate_symbol):
                        if value not in conflicts:
                            conflicts.append(value)
            assertions = current.local.setdefault("assertions", []) if current.local is not None else []
            if record not in assertions:
                assertions.append(record)
            for key in keys:
                key_to_indices.setdefault(key, set()).add(fused.index(current))
            continue
        if existing_indices:
            local_identifier_conflicts += 1
        item = FusedTarget(
            target_id=target_id,
            target_name=target_name or target_id,
            source="local",
            evidence_status="known_only",
            local={**record, "assertions": list(record.get("assertions", []))},
            alignment_status="local_only",
        )
        index = len(fused)
        fused.append(item)
        for key in keys:
            key_to_indices.setdefault(key, set()).add(index)

    excluded_non_human = 0
    excluded_unknown_species = 0
    excluded_species_conflict = 0
    excluded_non_protein = 0
    excluded_unknown_target_type = 0
    excluded_unreliable = 0
    ambiguous_alignments = 0
    duplicate_plato_records = 0
    accepted_plato = 0
    predicted_only_added = 0
    local_alignments_added = 0
    ordered_plato = sorted(plato_targets, key=_deterministic_target_key)
    for target in ordered_plato:
        if human_only and not target.is_human():
            excluded_non_human += 1
            if target.human_status == "unknown_species":
                excluded_unknown_species += 1
            elif target.human_status == "species_conflict":
                excluded_species_conflict += 1
            continue
        if not target.is_protein_target():
            excluded_non_protein += 1
            if not target.target_type:
                excluded_unknown_target_type += 1
            continue
        if reliable_only and target.reliable is not True:
            excluded_unreliable += 1
            continue

        matches = {
            index for key in target.alignment_keys for index in key_to_indices.get(key, set())
        }
        conflicting_match = len(matches) == 1 and _fused_identifier_conflict(fused[next(iter(matches))], target)
        if len(matches) > 1 or conflicting_match:
            ambiguous_alignments += 1
            if predicted_only_added >= max_plato_targets:
                continue
            item = FusedTarget(
                target_id=target.target_id,
                target_name=target.target_name,
                source="PLATO",
                evidence_status="predicted_only",
                plato=target,
                alignment_status="ambiguous_not_merged",
            )
            fused.append(item)
            accepted_plato += 1
            predicted_only_added += 1
            continue
        if len(matches) == 1:
            matching_index = next(iter(matches))
            current = fused[matching_index]
            current_symbols = {
                _normalize_identifier(_text(_first(current.local or {}, ("target_symbol", "gene_symbol", "symbol")))),
                _normalize_identifier(current.plato.target_symbol if current.plato else ""),
            } - {""}
            candidate_symbol = _normalize_identifier(target.target_symbol)
            if candidate_symbol and current_symbols and candidate_symbol not in current_symbols:
                target_symbol_disagreements += 1
            if current.plato is not None:
                duplicate_plato_records += 1
            if current.plato is None:
                accepted_plato += 1
                if current.local is not None:
                    local_alignments_added += 1
            if current.plato is not None:
                current.plato = _merge_duplicate_plato_group((current.plato, target))
                if current.local is None:
                    current.target_id = current.plato.target_id
                    current.target_name = current.plato.target_name
            else:
                current.plato = target
                if current.local is None:
                    current.target_id = target.target_id
                    current.target_name = target.target_name
            if current.local is not None:
                current.source = "local+PLATO"
                current.evidence_status = "known_and_predicted"
                preferred_symbol = _text(_first(current.local, ("target_symbol", "gene_symbol", "symbol")))
                preferred_symbol = preferred_symbol or _text(current.plato.target_symbol if current.plato else "")
                if preferred_symbol:
                    current.target_id = preferred_symbol.upper()
            current.alignment_method = next(
                (
                    key.split(":", 1)[0]
                    for key in target.alignment_keys
                    if matching_index in key_to_indices.get(key, set())
                ),
                current.alignment_method,
            )
            current.alignment_status = "aligned"
            for key in target.alignment_keys:
                key_to_indices.setdefault(key, set()).add(matching_index)
            continue

        if predicted_only_added >= max_plato_targets:
            continue
        item = FusedTarget(
            target_id=target.target_id,
            target_name=target.target_name,
            source="PLATO",
            evidence_status="predicted_only",
            plato=target,
            alignment_status="plato_only",
        )
        index = len(fused)
        fused.append(item)
        accepted_plato += 1
        predicted_only_added += 1
        for key in target.alignment_keys:
            key_to_indices.setdefault(key, set()).add(index)

    stats = {
        "local_targets": len(local_targets),
        "excluded_local_non_human": excluded_local_non_human,
        "excluded_local_non_protein": excluded_local_non_protein,
        "local_target_type_unspecified": local_target_type_unspecified,
        "local_target_type_unknown": local_target_type_unknown,
        "local_origin_unspecified": local_origin_unspecified,
        "local_identifier_conflicts": local_identifier_conflicts,
        "plato_targets_received": len(plato_targets),
        "plato_targets_accepted": accepted_plato,
        "plato_local_alignments_added": local_alignments_added,
        "plato_predicted_only_added": predicted_only_added,
        "plato_predicted_only_quota": max_plato_targets,
        "duplicate_plato_records": duplicate_plato_records,
        "ambiguous_alignments": ambiguous_alignments,
        "target_symbol_disagreements": target_symbol_disagreements,
        "excluded_non_human": excluded_non_human,
        "excluded_unknown_species": excluded_unknown_species,
        "excluded_species_conflict": excluded_species_conflict,
        "excluded_non_protein": excluded_non_protein,
        "excluded_unknown_target_type": excluded_unknown_target_type,
        "excluded_unreliable": excluded_unreliable,
        "known_and_predicted": sum(item.evidence_status == "known_and_predicted" for item in fused),
        "known_only": sum(item.evidence_status == "known_only" for item in fused),
        "predicted_only": sum(item.evidence_status == "predicted_only" for item in fused),
        "plato_applicability_known": sum(
            item.plato is not None and item.plato.applicability_score is not None for item in fused
        ),
        "total_fused": len(fused),
    }
    return fused, stats


class TargetResolver:
    """Resolve local evidence and cached PLATO predictions through one API."""

    def __init__(
        self,
        database: str | Path | None = None,
        plato_cache: str | Path | None = None,
        *,
        human_only: bool = True,
        reliable_only: bool = False,
        max_plato_targets: int = 20,
        provider_release: str = "r37",
        target_mapping: str | Path | None = None,
        target_mapping_sha256: str = "",
        ignore_stereochemistry: bool = False,
    ):
        self.database = Path(database) if database else None
        self.cache = PlatoPredictionCache(plato_cache) if plato_cache else None
        if not isinstance(ignore_stereochemistry, bool):
            raise TypeError("ignore_stereochemistry must be true or false")
        self.ignore_stereochemistry = ignore_stereochemistry
        self.target_mapper = ChemblTargetMapper(
            target_mapping,
            expected_sha256=target_mapping_sha256,
            release=provider_release,
        ) if target_mapping else ChemblTargetMapper(None, release=provider_release)
        self.human_only = human_only
        self.reliable_only = reliable_only
        if int(max_plato_targets) < 0:
            raise ValueError("max_plato_targets must be non-negative")
        self.max_plato_targets = int(max_plato_targets)
        self.provider_release = provider_release
        self.local_map: Dict[str, List[Any]] = {}
        self.last_local_resolution: Dict[str, Any] = {"status": "not_resolved"}
        self._database_hash_cache: Dict[Tuple[Any, ...], Tuple[str, str, List[str]]] = {}

    def _database_artifact_hash(self) -> Tuple[str, str, List[str]]:
        """Hash a stable SQLite snapshot, including a non-empty WAL when present."""
        if self.database is None:
            return "", "database_not_configured", []
        main = self.database.resolve()
        wal = Path(f"{main}-wal")
        components = [main]
        if wal.exists() and wal.stat().st_size > 0:
            components.append(wal)

        def signature(paths: Sequence[Path]) -> Tuple[Any, ...]:
            values: List[Any] = []
            for path in paths:
                stat = path.stat()
                values.extend((str(path), stat.st_size, stat.st_mtime_ns))
            return tuple(values)

        try:
            before = signature(components)
            cached = self._database_hash_cache.get(before)
            if cached is not None:
                return cached
            digest = sha256()
            for component in components:
                digest.update(component.name.encode("utf-8"))
                digest.update(b"\0")
                with component.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
            after = signature(components)
        except OSError:
            return "", "database_artifact_unreadable", [str(path) for path in components]
        if before != after:
            return "", "database_artifact_changed_during_hash", [str(path) for path in components]
        status = "sqlite_snapshot_with_wal" if len(components) == 2 else "sqlite_static_snapshot"
        result = (digest.hexdigest(), status, [str(path) for path in components])
        self._database_hash_cache = {before: result}
        return result

    def _query_database_snapshot(
        self, query: str, parameters: Tuple[Any, ...]
    ) -> Tuple[List[Tuple[Any, ...]], int, Dict[str, str], str, str, List[str]]:
        """Bind queried rows and artifact hash to one stable read transaction."""
        with closing(sqlite3.connect(self.database)) as connection:
            connection.execute("BEGIN")
            data_version_before = int(connection.execute("PRAGMA data_version").fetchone()[0])
            rows = connection.execute(query, parameters).fetchall()
            user_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
            database_metadata = _sqlite_resource_metadata(connection)
            artifact_hash, hash_status, hash_components = self._database_artifact_hash()
            data_version_after = int(connection.execute("PRAGMA data_version").fetchone()[0])
            connection.rollback()
        if data_version_before != data_version_after:
            return (
                rows, user_version, database_metadata, "",
                "database_changed_during_query_and_hash", hash_components,
            )
        return rows, user_version, database_metadata, artifact_hash, hash_status, hash_components

    def database_artifact_info(self) -> Dict[str, Any]:
        artifact_hash, status, components = self._database_artifact_hash()
        return {"sha256": artifact_hash, "status": status, "components": components}

    @staticmethod
    def _drug_key(value: str) -> str:
        return _normalize_name(value)

    def register_local_targets(self, drug_name: str, targets: Sequence[Any]) -> None:
        key = self._drug_key(drug_name)
        current = self.local_map.setdefault(key, [])
        for target in targets:
            record = _local_record(target)
            record["input_provenance_status"] = "registered_in_memory_unverified"
            if record.get("target_id") and record not in current:
                current.append(record)

    def local_targets(self, drug_name: str, smiles: str = "") -> List[Any]:
        registered = self.local_map.get(self._drug_key(drug_name))
        if registered:
            self.last_local_resolution = {
                "status": "registered_targets",
                "source": "in_memory",
                "records": len(registered),
            }
            return list(registered)
        if not self.database:
            self.last_local_resolution = {"status": "disabled", "records": 0}
            return []
        if not self.database.exists():
            self.last_local_resolution = {
                "status": "database_missing",
                "database_path": str(self.database),
                "records": 0,
            }
            return []
        query = """
            SELECT d.drug_id, d.drug_name, d.smiles,
                   t.target_symbol, t.target_name, t.target_type, t.uniprot_id,
                   t.source, t.external_id, g.taxid,
                   dt.source, dt.activity_type, dt.pchembl_value,
                   dt.action, dt.known_action, dt.evidence, dt.rowid
            FROM drug d
            JOIN drug_target dt ON dt.drug_id = d.drug_id
            JOIN target t ON t.target_id = dt.target_id
            LEFT JOIN gene g ON upper(trim(g.gene_symbol)) = upper(trim(t.target_symbol))
            WHERE lower(trim(d.drug_name)) = lower(trim(?))
               OR (? <> '' AND trim(d.smiles) = trim(?))
        """
        try:
            (
                rows, user_version, database_metadata, database_sha256,
                database_hash_status, database_hash_components,
            ) = self._query_database_snapshot(query, (drug_name, smiles, smiles))
            rich_schema = True
            rich_error = ""
        except sqlite3.Error as exc:
            rich_error = str(exc)
            missing_schema_element = "no such column" in rich_error.lower() or "no such table: gene" in rich_error.lower()
            if not missing_schema_element:
                self.last_local_resolution = {
                    "status": "database_error",
                    "database_path": str(self.database),
                    "error": rich_error,
                    "records": 0,
                }
                return []
            legacy_query = """
                SELECT d.drug_id, d.drug_name, d.smiles,
                       t.target_symbol, t.target_name, t.uniprot_id, t.external_id
                FROM drug d
                JOIN drug_target dt ON dt.drug_id = d.drug_id
                JOIN target t ON t.target_id = dt.target_id
                WHERE lower(trim(d.drug_name)) = lower(trim(?))
                   OR (? <> '' AND trim(d.smiles) = trim(?))
            """
            try:
                (
                    rows, user_version, database_metadata, database_sha256,
                    database_hash_status, database_hash_components,
                ) = self._query_database_snapshot(legacy_query, (drug_name, smiles, smiles))
                rich_schema = False
            except sqlite3.Error as legacy_exc:
                self.last_local_resolution = {
                    "status": "database_error",
                    "database_path": str(self.database),
                    "error": str(legacy_exc),
                    "rich_query_error": rich_error,
                    "records": 0,
                }
                return []
        name_match_ids = {
            row[0] for row in rows if _normalize_name(_text(row[1])) == _normalize_name(drug_name)
        }
        smiles_match_ids = {
            row[0] for row in rows if _text(smiles) and _text(row[2]) == _text(smiles)
        }
        matched_ids = name_match_ids | smiles_match_ids
        identity_conflict = bool(
            (_text(smiles) and name_match_ids and name_match_ids != smiles_match_ids)
            or len(matched_ids) > 1
        )
        if identity_conflict:
            self.last_local_resolution = {
                "status": "drug_identity_conflict",
                "database_path": str(self.database),
                "schema": "rich" if rich_schema else "legacy_fallback",
                "schema_user_version": user_version,
                "name_matched_drug_ids": sorted(name_match_ids),
                "smiles_matched_drug_ids": sorted(smiles_match_ids),
                "records": 0,
            }
            return []
        selected_drug_id = next(iter(matched_ids), None)
        if selected_drug_id is not None:
            rows = [row for row in rows if row[0] == selected_drug_id]
        grouped: Dict[str, Dict[str, Any]] = {}
        for row in rows:
            if rich_schema:
                _, _, _, symbol, name, target_type, uniprot, target_source, external_id, taxid, assertion_source, activity_type, pchembl, action, known_action, evidence, assertion_rowid = row
            else:
                _, _, _, symbol, name, uniprot, external_id = row
                target_type = target_source = taxid = assertion_source = activity_type = pchembl = action = known_action = evidence = ""
                assertion_rowid = None
            target_id = _text(symbol or uniprot or external_id or name).upper()
            if not target_id:
                continue
            record = grouped.setdefault(target_id, {
                "target_id": target_id,
                "target_name": _text(name or symbol or uniprot or external_id),
                "target_symbol": _text(symbol),
                "uniprot_id": _text(uniprot),
                "target_chembl_id": _text(external_id) if _text(external_id).upper().startswith("CHEMBL") else "",
                "external_id": _text(external_id),
                "target_type": _text(target_type),
                "taxid": "",
                "observed_taxids": [],
                "source": _text(target_source) or "local_database",
                "organism_status": "local_database_human_assumption",
                "input_provenance_status": "sqlite_artifact_audited" if rich_schema and database_sha256 else "sqlite_artifact_unverified",
                "provenance": {
                    "artifact_type": "sqlite",
                    "database_path": str(self.database.resolve()),
                    "database_sha256": database_sha256,
                    "hash_status": database_hash_status,
                    "hash_components": database_hash_components,
                    "schema_user_version": user_version,
                    "dataset_version": database_metadata.get("dataset_version", ""),
                    "mapping_version": database_metadata.get("mapping_version", ""),
                },
                "assertions": [],
            })
            normalized_taxid = _taxid(taxid)
            if normalized_taxid and normalized_taxid not in record["observed_taxids"]:
                record["observed_taxids"].append(normalized_taxid)
            assertion = {
                "source": _text(assertion_source),
                "activity_type": _text(activity_type),
                "pchembl_value": _float(pchembl),
                "action": _text(action),
                "known_action": _text(known_action),
                "evidence": _text(evidence),
                "stable_record_id": f"sqlite:drug_target:{assertion_rowid}" if assertion_rowid is not None else "",
                "provenance_hash": database_sha256,
            }
            if assertion not in record["assertions"]:
                record["assertions"].append(assertion)
        result = list(grouped.values())
        for record in result:
            taxids = sorted(set(record.get("observed_taxids", [])))
            record["observed_taxids"] = taxids
            if len(taxids) == 1:
                record["taxid"] = taxids[0]
                record["organism_status"] = "human_explicit" if taxids[0] == "9606" else "non_human"
            elif len(taxids) > 1:
                record["taxid"] = ""
                record["organism_status"] = "species_conflict"
        self.last_local_resolution = {
            "status": "ok" if result else "no_match",
            "database_path": str(self.database),
            "schema": "rich" if rich_schema else "legacy_fallback",
            "schema_user_version": user_version,
            "database_sha256": database_sha256,
            "database_hash_status": database_hash_status,
            "database_hash_components": database_hash_components,
            "dataset_version": database_metadata.get("dataset_version", ""),
            "mapping_version": database_metadata.get("mapping_version", ""),
            "rich_query_error": rich_error,
            "matched_by": (
                "drug_name_and_exact_smiles" if name_match_ids and smiles_match_ids
                else "drug_name" if name_match_ids else "exact_smiles" if smiles_match_ids else "none"
            ),
            "matched_drug_ids": sorted(matched_ids),
            "records": len(result),
        }
        return result

    def resolve(
        self,
        drug_name: str,
        smiles: str = "",
        local_targets: Optional[Sequence[Any]] = None,
    ) -> Tuple[List[FusedTarget], Dict[str, Any]]:
        if local_targets is not None:
            local = [
                {**_local_record(item), "input_provenance_status": "caller_supplied_unverified"}
                for item in local_targets
            ]
            local_resolution = {
                "status": "caller_supplied",
                "source": "explain_combination_argument",
                "records": len(local),
            }
        else:
            local = self.local_targets(drug_name, smiles)
            local_resolution = dict(self.last_local_resolution)
            if local_resolution.get("status") == "registered_targets":
                local = [
                    {**_local_record(item), "input_provenance_status": "registered_in_memory_unverified"}
                    for item in local
                ]
        plato: List[PlatoTargetEvidence] = []
        cache_meta: Dict[str, Any] = {
            "status": "disabled" if self.cache is None else ("cache_miss" if _text(smiles) else "not_applicable")
        }
        if self.cache is not None and _text(smiles):
            try:
                plato, cache_meta = (self.cache.load_ignoring_stereochemistry(smiles)
                                    if self.ignore_stereochemistry else self.cache.load(smiles))
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
                cache_meta = {"status": "cache_error", "error": str(exc)}
        target_mapping_meta: Dict[str, Any] = {"status": "disabled"}
        if plato:
            plato, target_mapping_meta = self.target_mapper.enrich(plato)
        elif self.target_mapper.status != "disabled":
            target_mapping_meta = {
                "status": self.target_mapper.status,
                "path": str(self.target_mapper.path) if self.target_mapper.path else "",
                "sha256": self.target_mapper.sha256,
                "mapped_targets": 0,
                "warnings": list(self.target_mapper.load_warnings),
            }
        fused, stats = fuse_local_and_plato_targets(
            local,
            plato,
            human_only=self.human_only,
            reliable_only=self.reliable_only,
            max_plato_targets=self.max_plato_targets,
        )
        actual_release = _text(cache_meta.get("provider_release")) or "unknown"
        release_mismatch = (
            actual_release != "unknown" and _text(self.provider_release)
            and actual_release.lower() != _text(self.provider_release).lower()
        )
        return fused, {
            "provider": "PLATO",
            "provider_release": actual_release,
            "expected_provider_release": self.provider_release,
            "release_mismatch": bool(release_mismatch),
            "cache": cache_meta,
            "target_mapping": target_mapping_meta,
            "fusion": stats,
            "local": local_resolution,
            "score_semantics": "PLATO score retained for ranking provenance; not treated as a probability",
        }
