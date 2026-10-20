"""Prepare a source review packet, without granting schema or raw-hash trust."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import sys
from urllib.request import Request, urlopen

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO.engineering import read_json, file_sha, selection_cases, verified_json
from InterpreterPLATO.research import stable_sha256, write_json_with_sha256
from InterpreterPLATO.targets import inspect_plato_payload, parse_plato_json, standardize_smiles

SITE_URLS = (
    "https://prometheus.farmacia.uniba.it/plato/pred.html",
    "https://prometheus.farmacia.uniba.it/plato/release.html",
    "https://prometheus.farmacia.uniba.it/plato/data.html",
)


def build_packet(selection_path, identity_path, raw_dir, output, capture_site=False):
    selection_path, identity_path, raw_dir, output = map(Path, (selection_path, identity_path, raw_dir, output))
    if output.exists():
        raise FileExistsError("Use a fresh packet directory; previous packets are retained")
    selection, identity = read_json(selection_path), verified_json(identity_path)
    if identity["selection_payload_sha256"] != stable_sha256(selection):
        raise ValueError("Identity audit is for a different selection")
    cases = selection_cases(selection)
    packet = {"contract": "plato_source_review_packet_v1", "created_at_utc": datetime.now(timezone.utc).isoformat(),
              "selection_id": selection["selection_id"], "selection_sha256": file_sha(selection_path),
              "identity_audit_path": str(identity_path.resolve()), "identity_audit_sha256": file_sha(identity_path),
              "review_status": "awaiting_source_and_schema_approval", "production_approval_granted": False,
              "trusted_raw_hashes_installed": False, "site_snapshots": [], "drugs": {}, "cases": []}
    output.mkdir(parents=True)
    if capture_site:
        for url in SITE_URLS:
            item = {"url": url}
            try:
                with urlopen(Request(url, headers={"User-Agent": "InterpreterPLATO/5.5 source-review"}), timeout=25) as response:
                    raw = response.read(2 * 1024 * 1024 + 1)
                    if len(raw) > 2 * 1024 * 1024:
                        raise ValueError("Oversized site snapshot")
                    if response.geturl() not in SITE_URLS:
                        raise ValueError("Unexpected final URL")
                    item.update(http_status=response.status, final_url=response.geturl())
                path = output / url.rsplit("/", 1)[-1]
                with path.open("xb") as handle:
                    handle.write(raw)
                item.update(path=str(path.resolve()), sha256=sha256(raw).hexdigest(), status="archived")
            except Exception as exc:
                item.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            packet["site_snapshots"].append(item)
    captures = {}
    for metadata_path in sorted(raw_dir.glob("*.metadata.json")):
        metadata = read_json(metadata_path)
        name = metadata["drug_name"]
        if name not in identity["drugs"]:
            continue
        if name in captures:
            raise ValueError(f"Multiple captures for {name}; explicitly select a raw input bundle")
        path = metadata_path.with_name(metadata_path.name.replace(".metadata.json", ".json"))
        if file_sha(path) != metadata["response_sha256"]:
            raise ValueError(f"Raw bytes changed for {name}")
        payload = read_json(path)
        parsed = parse_plato_json(payload)
        report = inspect_plato_payload(payload, parsed)
        capture_structure = standardize_smiles(metadata["submitted_smiles"])
        if capture_structure.get("inchikey") != identity["drugs"][name]["submitted_structure"].get("inchikey"):
            raise ValueError(f"Submitted query identity mismatch: {name}")
        captures[name] = {"raw_path": str(path.resolve()), "raw_sha256": file_sha(path),
                          "metadata_path": str(metadata_path.resolve()), "metadata_sha256": file_sha(metadata_path),
                          "captured_at_utc": metadata["captured_at_utc"], "endpoint": metadata["endpoint"],
                          "http_status": metadata["http_status"], "query_structure": capture_structure,
                          "release_from_payload": report["reported_release"],
                          "release_declared_by_capture": metadata["release"],
                          "query_and_release_not_embedded_in_raw_response": True,
                          "schema_fingerprint": report["schema_fingerprint"],
                          "schema_fixture_verified": report["schema_fixture_verified"],
                          "counts": {"parsed": len(parsed), "reliable_all_species": sum(p.reliable is True for p in parsed),
                                     "human": sum(p.is_human() for p in parsed),
                                     "reliable_human_protein_candidates_before_gate": sum(
                                         p.is_human() and p.is_protein_target() and p.reliable is True for p in parsed)},
                          "parser_audit": report,
                          "production_approval_granted": False}
    for name, record in identity["drugs"].items():
        packet["drugs"][name] = {"identity_status": record["status"], "identity_references": record.get("records", []),
                                  "unassigned_atom_stereocenters": record.get("unassigned_atom_stereocenters"),
                                  "capture": captures.get(name), "capture_status": "archived_pending_approval" if name in captures
                                  else "not_queried_identity_unresolved"}
    for case in cases:
        names = [case["input"][side + "_name"] for side in ("drug1", "drug2")]
        issues = ["schema_and_independent_raw_hash_approval_pending"]
        for name in names:
            item = packet["drugs"][name]
            if item["identity_status"] != "exact_inchikey_match":
                issues.append(f"chemical_identity_unresolved:{name}")
            capture = item["capture"]
            if not capture:
                issues.append(f"raw_response_missing:{name}")
            elif capture["counts"]["reliable_human_protein_candidates_before_gate"] == 0:
                issues.append(f"no_reliable_human_protein_targets:{name}")
        packet["cases"].append({"case_id": case["case_id"], "drugs": names, "cell_line": case["input"]["cell_line"],
                                "state": "blocked", "blocking_reasons": issues,
                                "prospective_after_approval": "still_blocked" if len(issues) > 1 else "requires_runtime_recheck"})
    packet["counts"] = {"selected_cases": len(cases), "unique_drugs": len(identity["drugs"]),
                        "archived_plato_captures": len(captures),
                        "identity_statuses": dict(Counter(r["status"] for r in identity["drugs"].values())),
                        "unique_candidate_schema_fingerprints": len({r["schema_fingerprint"] for r in captures.values()})}
    write_json_with_sha256(output / "source_review_packet.json", packet)
    return packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--identity-audit", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--capture-site", action="store_true", help="Archive three current official HTML pages by GET")
    args = parser.parse_args()
    packet = build_packet(args.selection, args.identity_audit, args.raw_dir, args.output_dir, args.capture_site)
    print(packet["counts"])


if __name__ == "__main__":
    main()
