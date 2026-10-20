"""Offline engineering verification of frozen ScoreRI development candidates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PACKAGE = Path(__file__).resolve().parents[1]
if str(PACKAGE.parent) not in sys.path:
    sys.path.insert(0, str(PACKAGE.parent))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", help="Prepare PENDING identity review records; never approve inputs")
    prepare.add_argument("--selection", required=True, type=Path)
    prepare.add_argument("--output", required=True, type=Path)
    for name in ("preflight", "run"):
        command = commands.add_parser(name)
        command.add_argument("--selection", required=True, type=Path)
        command.add_argument("--reviews", type=Path, help="Optional identity notes in test mode")
        command.add_argument("--input-policy", choices=("test", "strict"), default="test",
                             help="Default: original test flow without manual source acceptance")
        command.add_argument("--respect-stereochemistry", action="store_true",
                             help="Opt into exact stereochemical matching; default ignores stereo")
        command.add_argument("--context-config", type=Path, default=PACKAGE / "data/production_context.json")
        command.add_argument("--plato-cache", type=Path, default=PACKAGE / "data/processed/plato_cache")
        command.add_argument("--trusted-plato-sha256", action="append", default=[])
        command.add_argument("--output-dir", type=Path, required=True)
        if name == "run":
            command.add_argument("--variants", nargs="+", default=["full"])
            command.add_argument("--resume", action="store_true")
    summary = commands.add_parser("summarize")
    summary.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        from InterpreterPLATO.engineering import (
            preflight, read_json, review_template, run_suite, summarize,
        )
        from InterpreterPLATO.research import write_json_with_sha256
    except ImportError as exc:
        print(f"Environment dependency unavailable: {exc}. Use InterpreterPLATO/.venv/Scripts/python.exe (Windows), or the locked project environment.", file=sys.stderr)
        return 2
    try:
        if args.command == "prepare":
            if args.output.exists():
                raise ValueError("Review file exists; will not overwrite manual review records")
            write_json_with_sha256(args.output, review_template(read_json(args.selection)))
            print(f"Pending review template: {args.output}")
            return 0
        if args.command == "summarize":
            result = summarize(args.output_dir)
        else:
            pipeline, cases, report = preflight(args.selection, args.reviews, args.context_config,
                                               args.plato_cache, args.trusted_plato_sha256, args.input_policy,
                                               not args.respect_stereochemistry)
            if args.command == "preflight":
                destination = args.output_dir / "preflight.json"
                if destination.exists():
                    raise ValueError("Preflight report exists; use a new output directory to retain audit history")
                write_json_with_sha256(destination, report)
                print(f"Ready: {sum(case['ready'] for case in report['cases'])}/{len(cases)}; report: {destination}")
                for item in report["global_issues"]:
                    print(f"{item['code']}: {item['detail']}")
                return 0 if report["all_ready"] else 2
            print(f"Input policy: {args.input_policy}; ready: {sum(case['ready'] for case in report['cases'])}/{len(cases)}", flush=True)
            result = run_suite(args.output_dir, pipeline, cases, report, args.variants, args.resume)
        print(json.dumps({key: result[key] for key in ("complete", "n_planned", "counts", "problems")}, ensure_ascii=False))
        print(f"Verified summary: {args.output_dir / 'summary.json'}")
        return 0 if result["complete"] else 2
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
        print(f"Engineering verification stopped: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
