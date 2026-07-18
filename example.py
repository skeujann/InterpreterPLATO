"""Minimal production/default PLATO-only example.

The example reads manually archived PLATO responses; it does not query PLATO
over HTTP.  Until an r37 schema fixture is independently controlled and
verified, the production gate intentionally refuses to produce a score.
"""

from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from InterpreterPLATO import ExplanationPipeline


def main() -> int:
    package = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run one verified PLATO-only explanation")
    parser.add_argument("--drug1", required=True)
    parser.add_argument("--drug2", required=True)
    parser.add_argument("--drug1-smiles", required=True)
    parser.add_argument("--drug2-smiles", required=True)
    parser.add_argument("--cell-line", default="A549")
    parser.add_argument(
        "--context-config",
        default=str(package / "data" / "production_context.json"),
    )
    parser.add_argument(
        "--plato-cache",
        default=str(package / "data" / "processed" / "plato_cache"),
    )
    parser.add_argument("--plato-release", default="r37")
    parser.add_argument(
        "--trusted-plato-sha256",
        action="append",
        default=[],
        help="Independent raw PLATO export SHA-256; repeat for both drugs.",
    )
    args = parser.parse_args()

    pipe = ExplanationPipeline(
        config={
            "target_mode": "plato_only",
            "plato": {
                "cache_dir": args.plato_cache,
                "provider_release": args.plato_release,
                "trusted_export_sha256": args.trusted_plato_sha256,
            },
        }
    )
    pipe.load_production_resources(args.context_config)
    result = pipe.explain_combination(
        drug1_name=args.drug1,
        drug2_name=args.drug2,
        drug1_smiles=args.drug1_smiles,
        drug2_smiles=args.drug2_smiles,
        cell_line=args.cell_line,
        disease_label="pan_cancer_tier1",
    )
    result.print_summary()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
