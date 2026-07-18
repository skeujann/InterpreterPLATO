"""Run one user-triggered PLATO r37 target-fishing query and archive it."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any, Dict
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO.targets import PlatoPredictionCache, parse_plato_json, standardize_smiles


DEFAULT_ENDPOINT = "https://prometheus.farmacia.uniba.it/plato/pred.html"
ALLOWED_ENDPOINTS = {DEFAULT_ENDPOINT}
DEFAULT_MAX_RESPONSE_BYTES = 10 * 1024 * 1024
MAX_CONFIGURED_RESPONSE_BYTES = 64 * 1024 * 1024


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()
    return cleaned or "compound"


def _loads_strict_json(data: bytes) -> Any:
    return json.loads(
        data.decode("utf-8-sig"),
        parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError(f"non-finite JSON constant: {value}")
        ),
    )


def _validate_plato_endpoint(value: str, *, label: str = "endpoint") -> str:
    """Return the canonical endpoint or reject non-HTTPS/non-PLATO URLs."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty HTTPS URL")
    candidate = value.strip()
    try:
        parsed = urlsplit(candidate)
        explicit_port = parsed.port
    except ValueError as exc:
        raise ValueError(f"{label} is malformed: {exc}") from exc
    if (
        parsed.scheme.lower() != "https"
        or (parsed.hostname or "").lower() != "prometheus.farmacia.uniba.it"
        or explicit_port is not None
        or parsed.username is not None
        or parsed.password is not None
        or parsed.path != "/plato/pred.html"
        or parsed.query
        or parsed.fragment
    ):
        raise ValueError(
            f"{label} must be the allowlisted PLATO HTTPS target-fishing endpoint"
        )
    return DEFAULT_ENDPOINT


class _AllowlistedHTTPSRedirectHandler(HTTPRedirectHandler):
    """Reject a redirect before urllib contacts a non-allowlisted destination."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ANN001
        _validate_plato_endpoint(newurl, label="redirect target URL")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open_plato_request(request: Request, timeout: int):
    opener = build_opener(_AllowlistedHTTPSRedirectHandler())
    return opener.open(request, timeout=timeout)


def _read_response_limited(response: Any, max_response_bytes: int) -> bytes:
    content_length = response.headers.get("Content-Length", "")
    if content_length not in (None, ""):
        try:
            declared_length = int(str(content_length).strip())
        except (TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError("PLATO returned an invalid Content-Length header") from exc
        if declared_length < 0:
            raise RuntimeError("PLATO returned a negative Content-Length header")
        if declared_length > max_response_bytes:
            raise RuntimeError(
                "PLATO response exceeds the configured size limit "
                f"({declared_length} > {max_response_bytes} bytes)"
            )
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = response.read(min(64 * 1024, max_response_bytes - total + 1))
        if not chunk:
            break
        if not isinstance(chunk, bytes):
            raise RuntimeError("PLATO response stream returned non-byte data")
        total += len(chunk)
        if total > max_response_bytes:
            raise RuntimeError(
                f"PLATO response exceeded {max_response_bytes} bytes while streaming"
            )
        chunks.append(chunk)
    return b"".join(chunks)


def capture_plato_single_query(
    *,
    drug_name: str,
    smiles: str,
    endpoint: str = DEFAULT_ENDPOINT,
    release: str = "r37",
    timeout: int = 60,
    max_response_bytes: int = DEFAULT_MAX_RESPONSE_BYTES,
    raw_dir: str | Path | None = None,
    cache_dir: str | Path | None = None,
) -> Dict[str, Any]:
    """Submit one PLATO target-fishing form request and archive the response."""
    package = Path(__file__).resolve().parents[1]
    endpoint = _validate_plato_endpoint(endpoint)
    structure = standardize_smiles(smiles)
    if structure.get("valid") is not True:
        raise ValueError(
            "SMILES must be RDKit-validated before submitting to PLATO: "
            + ", ".join(structure.get("warnings", []))
        )
    if timeout < 5 or timeout > 300:
        raise ValueError("timeout must be between 5 and 300 seconds")
    if (
        isinstance(max_response_bytes, bool)
        or not isinstance(max_response_bytes, int)
        or max_response_bytes < 1
        or max_response_bytes > MAX_CONFIGURED_RESPONSE_BYTES
    ):
        raise ValueError(
            f"max_response_bytes must be an integer between 1 and {MAX_CONFIGURED_RESPONSE_BYTES}"
        )

    form = {
        "smiles": smiles,
        "tipoout": "json",
        "tipopred": "mussel",
        "invia": "get a response",
    }
    encoded = urlencode(form).encode("utf-8")
    request = Request(
        endpoint,
        data=encoded,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "InterpreterPLATO/5.1 single-query user-triggered",
        },
        method="POST",
    )
    try:
        with _open_plato_request(request, timeout) as response:
            status = int(getattr(response, "status", 0) or response.getcode())
            geturl = getattr(response, "geturl", None)
            observed_final_url = geturl() if callable(geturl) else endpoint
            final_endpoint = _validate_plato_endpoint(
                observed_final_url, label="redirected final URL",
            )
            if status != 200:
                raise RuntimeError(f"PLATO returned HTTP status {status}")
            raw = _read_response_limited(response, max_response_bytes)
            content_type = response.headers.get("Content-Type", "")
    except Exception as exc:  # noqa: BLE001 - external request failure is user-facing.
        raise RuntimeError(f"PLATO single query failed: {exc}") from exc

    try:
        payload = _loads_strict_json(raw)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"PLATO response was not strict JSON: {exc}") from exc
    parsed = parse_plato_json(payload)
    if not parsed:
        raise ValueError("PLATO response did not contain recognized target-fishing records")

    raw_sha256 = sha256(raw).hexdigest()
    captured_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    raw_root = Path(
        raw_dir
        if raw_dir is not None
        else package / "data" / "raw" / "plato" / "r37" / "single_query"
    ).expanduser()
    raw_root.mkdir(parents=True, exist_ok=True)
    base_name = f"{captured_at.replace(':', '').replace('+', 'Z')}_{_slug(drug_name)}_target_fishing"
    raw_path = raw_root / f"{base_name}.json"
    sidecar_path = raw_root / f"{base_name}.metadata.json"
    if raw_path.exists():
        raise FileExistsError(f"raw output already exists: {raw_path}")
    raw_path.write_bytes(raw)

    metadata: Dict[str, Any] = {
        "schema_version": "1.0",
        "provider": "PLATO",
        "release": release,
        "acquisition_method": "user_triggered_single_form_post",
        "endpoint": final_endpoint,
        "requested_endpoint": endpoint,
        "response_final_url": final_endpoint,
        "final_url_changed_from_request": (
            str(observed_final_url).strip() != endpoint
        ),
        "transport_security": "allowlisted_https_endpoint_with_tls_validation",
        "provider_content_authentication": "not_cryptographically_signed_by_provider",
        "http_method": "POST",
        "http_status": status,
        "content_type": content_type,
        "captured_at_utc": captured_at,
        "drug_name": drug_name,
        "submitted_smiles": smiles,
        "canonical_isomeric_smiles": structure.get("canonical_isomeric_smiles", ""),
        "inchikey": structure.get("inchikey", ""),
        "form": {
            "smiles": smiles,
            "tipoout": "json",
            "tipopred": "mussel",
        },
        "response_path": str(raw_path.resolve()),
        "response_sha256": raw_sha256,
        "response_size_bytes": len(raw),
        "response_size_limit_bytes": max_response_bytes,
        "parsed_targets": len(parsed),
        "reliable_targets": sum(target.reliable is True for target in parsed),
    }
    sidecar_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )

    cache_root = Path(
        cache_dir
        if cache_dir is not None
        else package / "data" / "processed" / "plato_cache"
    )
    cache = PlatoPredictionCache(cache_root)
    try:
        cache_path = cache.import_file(
            smiles=smiles,
            result_file=raw_path,
            drug_name=drug_name,
            provider_release=release,
            confirmed_real_export=True,
            query_binding_evidence=metadata,
        )
    except Exception as exc:  # noqa: BLE001 - report preserved raw evidence.
        raise RuntimeError(
            f"raw PLATO response saved at {raw_path}; metadata saved at {sidecar_path}; "
            f"cache import failed: {exc}"
        ) from exc
    return {
        "drug_name": drug_name,
        "smiles": smiles,
        "release": release,
        "raw_path": raw_path,
        "sidecar_path": sidecar_path,
        "cache_path": cache_path,
        "cache_dir": cache_root,
        "raw_sha256": raw_sha256,
        "final_endpoint": final_endpoint,
        "response_size_bytes": len(raw),
        "parsed_targets": len(parsed),
        "reliable_targets": metadata["reliable_targets"],
        "captured_at_utc": captured_at,
    }


def main() -> int:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Submit exactly one PLATO r37 target-fishing form request, save the raw JSON, "
            "and archive it in the PLATO cache with audited query binding."
        )
    )
    parser.add_argument("--drug-name", required=True)
    parser.add_argument("--smiles", required=True, help="Exact SMILES to submit to PLATO")
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    parser.add_argument("--release", default="r37")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument(
        "--max-response-bytes",
        type=int,
        default=DEFAULT_MAX_RESPONSE_BYTES,
        help="Hard streaming limit for the PLATO JSON response (default: 10 MiB)",
    )
    parser.add_argument(
        "--raw-dir",
        default=str(package / "data" / "raw" / "plato" / "r37" / "single_query"),
    )
    parser.add_argument(
        "--cache-dir",
        default=str(package / "data" / "processed" / "plato_cache"),
    )
    parser.add_argument("--print-report", action="store_true")
    args = parser.parse_args()

    try:
        result = capture_plato_single_query(
            drug_name=args.drug_name,
            smiles=args.smiles,
            endpoint=args.endpoint,
            release=args.release,
            timeout=args.timeout,
            max_response_bytes=args.max_response_bytes,
            raw_dir=args.raw_dir,
            cache_dir=args.cache_dir,
        )
    except Exception as exc:  # noqa: BLE001 - preserve raw capture, report cache issue.
        parser.error(str(exc))

    print(result["raw_path"])
    print(result["sidecar_path"])
    print(result["cache_path"])
    print(f"raw_sha256={result['raw_sha256']}")
    print(
        f"parsed_targets={result['parsed_targets']} "
        f"reliable_targets={result['reliable_targets']}"
    )
    if args.print_report:
        cache = PlatoPredictionCache(args.cache_dir)
        _, cache_metadata = cache.load(args.smiles)
        print(json.dumps(cache_metadata, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
