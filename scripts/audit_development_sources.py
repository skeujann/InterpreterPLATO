"""Acquire primary chemical identity records for a frozen development selection.

Read-only PubChem PUG REST GETs, sequential and rate-limited. Archives response
bytes and comparison findings; does NOT approve source reviews or PLATO schemas.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import time
from urllib.parse import quote
from urllib.request import Request, urlopen

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem
from InterpreterPLATO.engineering import read_json, selection_cases, file_sha
from InterpreterPLATO.research import stable_sha256, write_json_with_sha256
from InterpreterPLATO.targets import standardize_smiles


def compare_identity(smiles, properties):
    local = standardize_smiles(smiles)
    findings = {"submitted_structure": local, "records": []}
    for row in properties:
        reference_smiles = row.get("SMILES", row.get("IsomericSMILES", ""))
        reference = standardize_smiles(reference_smiles)
        key = row.get("InChIKey", "")
        if not reference.get("valid") or key != reference.get("inchikey"):
            status = "reference_record_inconsistent"
        elif not local.get("valid"):
            status = "invalid_input_structure"
        elif local.get("inchikey") == key:
            status = "exact_inchikey_match"
        elif local["inchikey"].split("-")[0] == key.split("-")[0]:
            status = "same_connectivity_different_stereo_or_isotope_layer"
        else:
            status = "different_connectivity_or_form"
        findings["records"].append({"cid": row.get("CID"), "reference_inchikey": key,
                                    "reference_smiles": reference_smiles, "reference_structure": reference,
                                    "status": status, "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{row.get('CID')}"})
    mol = Chem.MolFromSmiles(smiles)
    findings["unassigned_atom_stereocenters"] = (
        sum(label == "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True)) if mol else None)
    findings["status"] = (findings["records"][0]["status"] if len(findings["records"]) == 1
                          else "ambiguous_or_empty_name_lookup")
    findings["review_approval_granted"] = False
    return findings


def acquire(selection_path, output):
    selection_path, output = Path(selection_path).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Use a new output directory; raw identity records are immutable")
    selection = read_json(selection_path)
    cases = selection_cases(selection)
    drugs = {}
    for case in cases:
        for side in ("drug1", "drug2"):
            drugs[case["input"][side + "_name"]] = case["input"][side + "_smiles"]
    output.mkdir(parents=True)
    report = {"contract": "chemical_source_comparison_v1", "selection_id": selection["selection_id"],
              "selection_sha256": file_sha(selection_path), "selection_payload_sha256": stable_sha256(selection),
              "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
              "source": "PubChem PUG REST", "no_plato_query_performed": True,
              "approval_granted": False, "drugs": {}}
    for index, (name, smiles) in enumerate(drugs.items(), 1):
        # This spelling maps to the same registered compound name, not an
        # inferred stereochemical replacement of the frozen input structure.
        query_name = "5-azacytidine" if name == "Azacytidine, 5-" else name
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + quote(query_name, safe="")
               + "/property/InChIKey,IsomericSMILES/JSON")
        record = {"dataset_id": name, "query_name": query_name, "query_url": url}
        try:
            request = Request(url, headers={"User-Agent": "InterpreterPLATO/5.5 chemical-source-audit"})
            with urlopen(request, timeout=25) as response:
                raw = response.read(2 * 1024 * 1024 + 1)
                if len(raw) > 2 * 1024 * 1024:
                    raise ValueError("Oversized PubChem response")
                record.update(http_status=response.status, response_url=response.geturl())
            raw_path = output / f"drug_{index:02d}.pubchem.json"
            with raw_path.open("xb") as handle:
                handle.write(raw)
            record.update(raw_file=raw_path.name, raw_sha256=sha256(raw).hexdigest())
            payload = json.loads(raw)
            record.update(compare_identity(smiles, payload["PropertyTable"]["Properties"]))
        except Exception as exc:
            record.update(status="acquisition_or_parse_failed", error=f"{type(exc).__name__}: {exc}")
        report["drugs"][name] = record
        write_json_with_sha256(output / "identity_audit.json", report)
        print(f"{index}/{len(drugs)} {name}: {record['status']}", flush=True)
        if index < len(drugs):
            time.sleep(0.3)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    acquire(args.selection, args.output_dir)


if __name__ == "__main__":
    main()
