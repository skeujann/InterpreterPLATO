"""Render saved result(s) without querying PLATO or rerunning inference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Optional, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from InterpreterPLATO.reporting import companion_report_path, write_readable_report
else:
    from ..reporting import companion_report_path, write_readable_report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Saved explanation, list of explanations, or explicit failure JSON")
    parser.add_argument("--output", help="Chinese Markdown report (defaults to INPUT.report.md)")
    args = parser.parse_args(argv)
    source = Path(args.input).resolve()
    destination = Path(args.output).resolve() if args.output else companion_report_path(source)
    if destination == source or destination.suffix.lower() != ".md":
        parser.error("--output must be a separate .md report, never the input JSON")
    try:
        payload = json.loads(source.read_text(encoding="utf-8-sig"))
        record = write_readable_report(destination, payload)
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        parser.error(f"unable to render saved explanation: {exc}")
    print(record["report_path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
