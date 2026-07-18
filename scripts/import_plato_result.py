"""Archive one manually obtained PLATO JSON response in the audit cache."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO.targets import PlatoPredictionCache


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Archive one manually obtained PLATO target-fishing JSON. "
            "Importing or attesting a file does not verify it as a PLATO response."
        )
    )
    parser.add_argument("--smiles", required=True, help="The exact SMILES submitted to PLATO")
    parser.add_argument("--result", required=True, help="Downloaded PLATO JSON result")
    parser.add_argument("--drug-name", default="")
    parser.add_argument(
        "--release",
        default="unknown",
        help="Release stated by PLATO; use unknown if the downloaded payload does not identify it",
    )
    parser.add_argument(
        "--attest-downloaded-from-plato",
        action="store_true",
        help=(
            "Record the user's statement that the file was downloaded from PLATO. "
            "This is not provider authentication or schema verification; without it "
            "the file is retained only as an unattested archive."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path(__file__).resolve().parents[1] / "data" / "processed" / "plato_cache"),
    )
    parser.add_argument("--print-report", action="store_true", help="Print cache/parser provenance after import")
    args = parser.parse_args()

    source = Path(args.result).expanduser()
    if not source.is_file():
        parser.error(f"PLATO result is not a readable file: {source}")
    try:
        json.loads(
            source.read_text(encoding="utf-8-sig"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-finite JSON constant: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        parser.error(f"PLATO result is not strict UTF-8 JSON: {exc}")

    cache = PlatoPredictionCache(args.cache_dir)
    try:
        # `confirmed_real_export` is the cache API's legacy keyword.  The cache
        # records it only as an unauthenticated user attestation.
        destination = cache.import_file(
            smiles=args.smiles,
            result_file=source,
            drug_name=args.drug_name,
            provider_release=args.release,
            confirmed_real_export=args.attest_downloaded_from_plato,
        )
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        parser.error(f"unable to archive PLATO result: {exc}")
    print(destination)
    if args.attest_downloaded_from_plato:
        print(
            "USER ATTESTATION RECORDED: this is not PLATO authentication, "
            "schema verification, or permission to produce a score."
        )
    else:
        print(
            "ARCHIVE ONLY: no source attestation was recorded; this cache entry "
            "cannot satisfy the production gate."
        )
    if args.print_report:
        _, metadata = cache.load(args.smiles)
        print(json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
