"""Previously authorized real-input test policy, isolated from strict execution.

Capture hashes freeze the bytes used for a test; they do not grant source approval.
Only the six historical approval requirements may be waived during validation.
"""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path

from .interpreter import ExplanationPipeline
from .targets import inspect_plato_payload, _loads_strict_json, _single_query_binding_report


APPROVAL_FAILURES = frozenset({
    "PLATO_r37_schema_fixture_not_verified", "generic_or_missing_PLATO_adapter",
    "PLATO_provider_not_schema_verified", "PLATO_target_fishing_method_not_schema_verified",
    "PLATO_export_sha256_not_independently_trusted", "PLATO_release_not_payload_verified",
})
OBSERVED_SCHEMAS = frozenset({
    "749895127c8b896ffbef42d0cb693dc018262a3febb637155d59845c774c5aa6",
    "b7fdb3e7b40bf1d26397ceb0182bb122d2f04970797efaf9749b4a11921d3403",
})


def test_policy(cache_inventory):
    return {
        "policy_id": "development_test_without_manual_acceptance_v1",
        "run_mode": "test", "independent_review_required_for_this_run": False,
        "independent_review_completed": False, "production_schema_approved": False,
        "capture_integrity_pins": dict(cache_inventory),
        "hash_semantics": "frozen input bytes, not independent source approval",
        "waivable_failures": sorted(APPROVAL_FAILURES),
        "observed_schema_fingerprints": sorted(OBSERVED_SCHEMAS),
        "structure_scope": "configured structure matching; exact original POST provenance retained",
        "credibility_semantics": "original strict source status retained",
    }


class DevelopmentTestPipeline(ExplanationPipeline):
    def __init__(self, config):
        if config.get("run_mode") != "test":
            raise ValueError("Source-review exception requires test mode")
        super().__init__(config=config)
        self._test_validation_active = False
        self.review_waivers = {}

    def _plato_cache_integrity_failures(self, resolution, expected_release,
                                       trusted_export_sha256=(), expected_structure=None):
        failures = super()._plato_cache_integrity_failures(
            resolution, expected_release, trusted_export_sha256, expected_structure)
        # Credibility and source reporting always see the original strict status.
        if not self._test_validation_active:
            return failures
        if self.config.get("run_mode") != "test":
            return failures + ["test_mode_changed"]
        cache = resolution.get("cache", {})
        source = cache.get("source", {})
        raw_hash = source.get("sha256", "")
        root = Path(self.config["plato"]["cache_dir"]).resolve()
        raw_path = (root / source.get("raw_blob_relative_path", "")).resolve()
        pins = self.config.get("test_source_review_exception", {}).get("capture_integrity_pins", {})
        if root not in raw_path.parents or pins.get(str(raw_path)) != raw_hash or not raw_hash:
            return failures + ["test_raw_input_not_pinned_or_changed"]
        try:
            raw_bytes = raw_path.read_bytes()
            if sha256(raw_bytes).hexdigest() != raw_hash:
                return failures + ["test_raw_input_not_pinned_or_changed"]
            inspection = inspect_plato_payload(_loads_strict_json(raw_bytes.decode("utf-8-sig")))
        except (OSError, ValueError, UnicodeError):
            return failures + ["test_raw_input_unreadable"]
        if (expected_release != "r37"
                or inspection["schema_fingerprint"] not in OBSERVED_SCHEMAS
                or cache.get("schema_fingerprint") != inspection["schema_fingerprint"]
                or inspection["field_names"] != ["chemblid", "nome", "score", "simili", "tid", "validato"]
                or inspection["provider_validation_status"] not in ("missing", "valid_plato")
                or inspection["method_validation_status"] not in ("missing", "valid_target_fishing")):
            return failures + ["test_schema_or_release_changed"]
        query = cache.get("query", {})
        # Group matching never changes which stereospecific structure was sent.
        submitted = query.get("canonical_isomeric_smiles", "")
        binding = _single_query_binding_report(submitted, raw_hash, query.get("single_query_evidence"))
        if query.get("query_binding_status") != "plato_single_query_post_verified" or not binding:
            return failures + ["test_official_request_binding_missing"]
        remaining = [value for value in failures if value not in APPROVAL_FAILURES]
        self.review_waivers[raw_hash] = {
            "waived_for_test": sorted(set(failures) & APPROVAL_FAILURES),
            "remaining_failures": remaining,
        }
        return remaining

    def _validate_plato_only_resolution(self, *args, **kwargs):
        self._test_validation_active = True
        try:
            return super()._validate_plato_only_resolution(*args, **kwargs)
        finally:
            self._test_validation_active = False

    def explain_combination(self, *args, **kwargs):
        self.review_waivers = {}
        result = super().explain_combination(*args, **kwargs)
        result.metadata["test_source_review_exception"] = {
            **deepcopy(self.config["test_source_review_exception"]),
            "applied_waivers": deepcopy(self.review_waivers),
        }
        result.metadata["formal_result_eligibility_reasons"].append(
            "development_test_without_manual_source_acceptance")
        return result
