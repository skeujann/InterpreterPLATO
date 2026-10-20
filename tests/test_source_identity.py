from InterpreterPLATO.scripts.audit_development_sources import compare_identity
from InterpreterPLATO.targets import standardize_smiles
from InterpreterPLATO.research import stable_sha256, write_json_with_sha256
from InterpreterPLATO.engineering import file_sha
from InterpreterPLATO.scripts.build_plato_review_packet import build_packet
from InterpreterPLATO.scripts.test_support import plato_payload
import pytest


def row(smiles):
    return {"CID": 1, "SMILES": smiles, "InChIKey": standardize_smiles(smiles)["inchikey"]}


def test_exact_match_is_not_approval():
    result = compare_identity("OCC", [row("CCO")])
    assert result["status"] == "exact_inchikey_match"
    assert result["review_approval_granted"] is False


def test_missing_stereo_is_not_promoted_to_matching_reference():
    result = compare_identity("FC(Cl)Br", [row("F[C@H](Cl)Br")])
    assert result["status"] == "same_connectivity_different_stereo_or_isotope_layer"
    assert result["unassigned_atom_stereocenters"] == 1


def test_mismatch_inconsistent_and_ambiguous_records_fail():
    assert compare_identity("CCN", [row("CCO")])["status"] == "different_connectivity_or_form"
    assert compare_identity("CCO", [row("CCO"), row("CCO")])["status"] == "ambiguous_or_empty_name_lookup"
    inconsistent = {**row("CCN"), "SMILES": "CCO"}
    assert compare_identity("CCO", [inconsistent])["status"] == "reference_record_inconsistent"


def packet_inputs(tmp_path):
    selection = {"selection_id": "source-test", "drug_structures": {
        "raw_smiles_by_dataset_id": {"a": "CCO", "b": "CCN"}}, "selected": [{
            "candidate_id": "DEV-1", "drug_a": "a", "drug_b": "b", "cell_line_name": "TEST",
            "pred_syn_mean": 2., "role": "development_only"}]}
    path = tmp_path / "selection.json"
    write_json_with_sha256(path, selection)
    audit = {"selection_payload_sha256": stable_sha256(selection), "drugs": {
        name: compare_identity(smiles, [row(smiles)]) for name, smiles in (("a", "CCO"), ("b", "CCN"))}}
    audit_path = tmp_path / "identity.json"
    write_json_with_sha256(audit_path, audit)
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    raw_path = raw_dir / "a.json"
    write_json_with_sha256(raw_path, plato_payload("CCO"))
    write_json_with_sha256(raw_dir / "a.metadata.json", {
        "drug_name": "a", "submitted_smiles": "CCO", "response_sha256": file_sha(raw_path),
        "captured_at_utc": "2026-10-13T00:00:00+00:00", "endpoint": "synthetic-fixture-not-network",
        "http_status": 200, "release": "r37"})
    return path, audit_path, raw_dir


def test_packet_does_not_approve_and_preserves_missing_capture(tmp_path):
    args = packet_inputs(tmp_path)
    packet = build_packet(*args, tmp_path / "packet")
    assert packet["production_approval_granted"] is False
    assert packet["trusted_raw_hashes_installed"] is False
    assert packet["counts"]["archived_plato_captures"] == 1
    assert packet["cases"][0]["state"] == "blocked"
    assert "raw_response_missing:b" in packet["cases"][0]["blocking_reasons"]
    assert not packet["site_snapshots"]


def test_packet_rejects_raw_tampering(tmp_path):
    args = packet_inputs(tmp_path)
    (args[2] / "a.json").write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="Raw bytes changed"):
        build_packet(*args, tmp_path / "packet")


def test_packet_rejects_query_structure_mismatch(tmp_path):
    args = packet_inputs(tmp_path)
    metadata_path = args[2] / "a.metadata.json"
    import json
    metadata = json.loads(metadata_path.read_text())
    metadata["submitted_smiles"] = "CCC"
    write_json_with_sha256(metadata_path, metadata)
    with pytest.raises(ValueError, match="query identity mismatch"):
        build_packet(*args, tmp_path / "packet")
