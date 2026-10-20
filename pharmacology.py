"""Drug-specific literature roles for existing PLATO candidates.

This annotation layer does not add targets, reorder context scores, establish
isoform binding from a family claim, or change statistical support.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Dict, Mapping, Sequence

from .context import ContextScore, gene
from .targets import standardize_smiles

ROLE_POLICY_VERSION = "drug_specific_literature_roles_v1"
# Pin the LF bytes stored in Git, as enforced by .gitattributes on checkout.
DEFAULT_ROLE_SHA256 = "f5a42b7b2f7d55d8b2438d9d9c8514d9f53ed04a67f183e65cd74ca69971be5f"
HANDLING_ROLES = frozenset({"drug_metabolism", "drug_transport", "drug_activation_retention"})
ALLOWED_ROLES = HANDLING_ROLES | {"pharmacodynamic_target", "context_association"}


class PharmacologyRoleCatalog:
    """Load one pinned, validated annotation resource; no online inference."""

    def __init__(self, path: str | Path, expected_sha256: str):
        self.path = Path(path).resolve()
        if not re.fullmatch(r"[0-9a-f]{64}", str(expected_sha256)):
            raise ValueError("pharmacology_roles.sha256 must be an explicit SHA-256")
        raw = self.path.read_bytes()
        self.sha256 = sha256(raw).hexdigest()
        if self.sha256 != expected_sha256:
            raise ValueError("Pharmacology role resource SHA-256 mismatch")
        payload = json.loads(raw.decode("utf-8-sig"))
        if not isinstance(payload, dict) or payload.get("schema_version") != "1.0":
            raise ValueError("Unsupported pharmacology role schema")
        if not all(isinstance(payload.get(k), str) and payload[k].strip() for k in ("version", "curation_status", "coverage")):
            raise ValueError("Role resource requires version, curation_status and coverage")
        references = payload.get("references")
        records = payload.get("records")
        if not isinstance(references, dict) or not isinstance(records, list):
            raise ValueError("Role resource requires references and records")
        for ref in references.values():
            if not isinstance(ref, dict) or not ref.get("title") or not str(ref.get("url", "")).startswith("https://"):
                raise ValueError("Each literature reference requires title and HTTPS URL")
        self._index: Dict[tuple, list] = {}
        seen_ids = set()
        seen_relations = set()
        for record in records:
            if not isinstance(record, dict):
                raise ValueError("Role record must be an object")
            required = ("id", "drug_name", "drug_inchikey", "target_gene", "role", "evidence_scope", "action", "process", "claim", "experimental_context")
            if not all(isinstance(record.get(k), str) and record[k].strip() for k in required):
                raise ValueError("Role record is missing a required field")
            if record["role"] not in ALLOWED_ROLES or record["evidence_scope"] not in {"exact_protein", "target_family", "association_only"}:
                raise ValueError("Unknown pharmacology role or evidence scope")
            if record["evidence_scope"] == "association_only" and record["role"] != "context_association":
                raise ValueError("An association alone cannot establish a pharmacology role")
            if record["role"] == "context_association" and record["evidence_scope"] != "association_only":
                raise ValueError("Context association must use association_only scope")
            if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", record["drug_inchikey"]):
                raise ValueError("Role record requires a full structure-specific InChIKey")
            if not re.fullmatch(r"[A-Z][A-Z0-9-]*", record["target_gene"]):
                raise ValueError("Role target must be an exact HGNC-style symbol, not a wildcard")
            refs = record.get("reference_ids")
            if not isinstance(refs, list) or not refs or any(not isinstance(ref, str) or ref not in references for ref in refs):
                raise ValueError("Role record has missing or unknown references")
            relation = tuple(record[k] for k in ("drug_inchikey", "target_gene", "role", "evidence_scope", "action"))
            if record["id"] in seen_ids or relation in seen_relations:
                raise ValueError("Duplicate role record or relationship")
            seen_ids.add(record["id"])
            seen_relations.add(relation)
            stored = deepcopy(record)
            stored["references"] = [{"id": ref, **deepcopy(references[ref])} for ref in dict.fromkeys(refs)]
            self._index.setdefault((record["drug_inchikey"], record["target_gene"]), []).append(stored)
        self._audit = {
            "policy_version": ROLE_POLICY_VERSION,
            "resource_path": str(self.path),
            "resource_sha256": self.sha256,
            "resource_version": payload["version"],
            "curation_status": payload["curation_status"],
            "coverage": payload["coverage"],
            "record_count": len(records),
            "drug_structure_count": len({key[0] for key in self._index}),
            "identity_policy": "full_standardized_inchikey_and_exact_candidate_gene",
            "target_source": "PLATO_only",
            "annotation_source": "literature_resource",
            "used_for_scoring_or_statistical_selection": False,
            "independent_validation_performed": False,
        }

    def audit(self) -> Dict[str, Any]:
        return deepcopy(self._audit)

    def annotate(self, profile: Sequence[ContextScore], smiles: str) -> list:
        identity = standardize_smiles(smiles)
        key = str(identity.get("inchikey") or "")
        result = []
        for target in profile:
            records = deepcopy(self._index.get((key, gene(target.target_id)), []))
            roles = sorted({record["role"] for record in records})
            annotation = {
                "policy_version": ROLE_POLICY_VERSION,
                "status": "literature_annotated" if records else "unresolved_role",
                "roles": roles or ["unknown"],
                "drug_inchikey": key,
                "identity_match": "exact_structure_and_gene" if records else "no_annotation_match",
                "resource_sha256": self.sha256,
                "records": records,
                "current_cell_line_validated": False,
                "combination_mechanism_validated": False,
            }
            result.append(replace(target, pharmacology=annotation))
        return result


def summarize_roles(profile: Sequence[ContextScore]) -> Dict[str, Any]:
    groups = {key: [] for key in ("direct_pharmacodynamic", "family_pharmacodynamic", "drug_handling", "context_association", "unresolved")}
    for item in profile:
        annotation = item.pharmacology
        records = annotation.get("records", [])
        entry = {"target_id": item.target_id, "target_name": item.target_name,
                 "context_score": item.final_score, "roles": annotation.get("roles", ["unknown"]),
                 "record_ids": [record["id"] for record in records]}
        if not records:
            groups["unresolved"].append(entry)
        if any(r["role"] == "pharmacodynamic_target" and r["evidence_scope"] == "exact_protein" for r in records):
            groups["direct_pharmacodynamic"].append(entry)
        if any(r["role"] == "pharmacodynamic_target" and r["evidence_scope"] == "target_family" for r in records):
            groups["family_pharmacodynamic"].append(entry)
        if any(r["role"] in HANDLING_ROLES for r in records):
            groups["drug_handling"].append(entry)
        if any(r["role"] == "context_association" for r in records):
            groups["context_association"].append(entry)
    return {"groups": groups, "counts": {k: len(v) for k, v in groups.items()},
            "multiple_roles_per_candidate_allowed": True,
            "unknown_does_not_mean_no_pharmacological_effect": True}


def role_composition(genes: Sequence[str], profile: Sequence[ContextScore]) -> Dict[str, Any]:
    selected = {gene(value) for value in genes}
    relevant = [item for item in profile if gene(item.target_id) in selected]
    summary = summarize_roles(relevant)
    unknown = sorted(selected - {gene(item.target_id) for item in relevant})
    all_handling = bool(selected) and not unknown and all(
        set(item.pharmacology.get("roles", [])) and
        set(item.pharmacology["roles"]) <= HANDLING_ROLES for item in relevant
    )
    return {"counts": summary["counts"], "unresolved_identifiers": unknown,
            "drug_handling_only": all_handling,
            "pharmacodynamic_enrichment_tested": False,
            "interpretation": "drug_handling_annotation_only" if all_handling else "mixed_or_unresolved_candidate_context"}


def pathway_role_context(enrichment: Mapping[str, Any], left: Sequence[ContextScore], right: Sequence[ContextScore]) -> Dict[str, Any]:
    return {
        side: [{"pathway": row["pathway"], "q": row["q"],
                **role_composition(row.get("hits", []), profile)}
               for row in enrichment.get(key, [])]
        for side, key, profile in (("drug1", "significant_drug1", left), ("drug2", "significant_drug2", right))
    }


def pharmacology_overview(left: Sequence[ContextScore], right: Sequence[ContextScore], names: Sequence[str]) -> Dict[str, Any]:
    summaries = {"drug1": summarize_roles(left), "drug2": summarize_roles(right)}
    clauses = []
    for name, profile in zip(names, (left, right)):
        anchors = []
        for item in profile:
            for record in item.pharmacology.get("records", []):
                if record["role"] == "pharmacodynamic_target":
                    scope = "family evidence only; isoform engagement unverified" if record["evidence_scope"] == "target_family" else "protein-specific literature"
                    anchors.append(f"{item.target_id}: {record['action']} / {record['process']} ({scope})")
        clauses.append(f"{name}: " + ("; ".join(anchors) if anchors else "no literature-assigned pharmacodynamic candidate in the available profile"))
    annotated = any(item.pharmacology.get("records") for item in list(left) + list(right))
    narrative = (
        "Literature pharmacology among existing PLATO candidates. " + ". ".join(clauses)
        + ". Metabolism, transport and intracellular retention annotations are reported separately. "
          "These annotations do not establish the combination mechanism in the current cell line."
    ) if annotated else "No matched drug-specific literature roles; the candidate-context interpretation remains exploratory."
    return {"policy_version": ROLE_POLICY_VERSION, "status": "partial_literature_coverage" if annotated else "unresolved",
            **summaries, "narrative": narrative,
            "family_evidence_is_isoform_validation": False,
            "roles_change_context_ranking": False,
            "roles_change_statistical_target_sets": False,
            "roles_add_score_bonus": False,
            "combination_mechanism_validated": False}


def summary_target_label(profile: Sequence[ContextScore], limit: int = 3) -> str:
    groups = summarize_roles(profile)["groups"]
    labels = [r["target_id"] for r in groups["direct_pharmacodynamic"]]
    labels += [r["target_id"] + " [family evidence]" for r in groups["family_pharmacodynamic"]]
    return ", ".join(labels[:limit]) if labels else "no assigned pharmacodynamic target"
