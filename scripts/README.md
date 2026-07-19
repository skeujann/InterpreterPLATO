# Script Entry Points

- `build_context_artifacts.py`: build the human cancer context from pinned HGNC, STRING, DepMap, Reactome, and COSMIC source files.
- `build_chembl_target_mapping.py`: build the ChEMBL r37 target-ID to HGNC-symbol translation layer; it does not import ChEMBL drug-target evidence.
- `verify_manifest.py`: verify file size and SHA-256 entries in a manifest.
- `query_pair_and_explain.py`: sequentially submit one PLATO target-fishing request per drug in a pair, then run the strict PLATO-only mechanism interpreter.
- `query_plato_once.py`: submit exactly one user-triggered PLATO target-fishing form request, save raw JSON plus a metadata sidecar, and import it with audited query binding.
- `import_plato_result.py`: archive one manually obtained PLATO JSON response. It does not contact PLATO and does not authenticate the platform.
- `run_interpreter.py`: strict PLATO-only single-pair or batch CLI.
- `test_*.py`: offline regression and serialization checks.

## Pair-Level Query Then Explain

```powershell
python InterpreterPLATO/scripts/query_pair_and_explain.py `
  --drug1 "Drug-A" `
  --drug2 "Drug-B" `
  --smiles-lookup C:\datasets\drug_smiles.csv `
  --cell-line A549 `
  --disease-label pan_cancer_tier1 `
  --quick
```

When `--drug1-smiles` or `--drug2-smiles` is omitted, an explicit repeatable
`--smiles-lookup PATH` is required. No SMILES table is bundled, and
InterpreterPLATO never scans ambient sibling projects. Alternatively, provide
both SMILES explicitly. Test resolution without contacting PLATO:

```powershell
python InterpreterPLATO/scripts/query_pair_and_explain.py `
  --drug1 "Gemcitabine" `
  --drug2 "Paclitaxel" `
  --smiles-lookup C:\datasets\drug_smiles.csv `
  --resolve-only
```

If one name maps to multiple distinct SMILES, the script fails closed and asks
for the explicit SMILES.

The script submits two PLATO form requests sequentially, with a configurable
`--delay-seconds` pause between them. It saves run-scoped raw responses,
sidecars, cache entries, and `run_manifest.json`. Same-run hashes are integrity
checks, not trust anchors. `explanation.json` is written only when caller-
supplied, independently approved `--trusted-plato-sha256` values match both
captures and all other gates pass; otherwise the evidence is preserved and the
command exits with code 2. The PLATO page path is `pred.html`, but the submitted
form requests `tipoout=json`, so the archived payload is JSON.

## ChEMBL Target Mapping

After placing the ChEMBL r37 SQLite archive and companion files under
`InterpreterPLATO/data/raw/target_mapping/chembl/r37/`, build the translation
layer:

```powershell
python InterpreterPLATO/scripts/build_chembl_target_mapping.py
```

The output maps PLATO `target_chembl_id` values to HGNC symbols for
PPI/expression/pathway/disease context lookup. It is not a local target source
and does not use ChEMBL drug-target associations.

## Single PLATO Query

```powershell
python InterpreterPLATO/scripts/query_plato_once.py `
  --drug-name "Drug-A" `
  --smiles "exact SMILES submitted to PLATO" `
  --release r37 `
  --print-report
```

One invocation sends one `POST` request to the official PLATO target-fishing
form endpoint. The script is intentionally not a generic client: it does not
loop, does not run concurrently, and rejects non-PLATO endpoints. It records
the raw response, SHA-256, submitted form, endpoint, HTTP method, and SMILES
identity evidence so the cache can use
`query_binding_status=plato_single_query_post_verified`.

## Manual Archive Is Not Verification

```powershell
python InterpreterPLATO/scripts/import_plato_result.py `
  --drug-name "Drug-A" `
  --smiles "exact SMILES submitted to PLATO" `
  --result "C:/path/drug-a-plato.json" `
  --release r37 `
  --attest-downloaded-from-plato `
  --print-report
```

The attestation records a user statement only. It is not provider
authentication, schema verification, or permission to produce a score.
Without a frozen verified r37 fixture, an imported response remains an audit
archive and the production CLI must reject scoring.
The current production schema allowlist is empty, so this restriction applies to
every imported response; the scoring command below documents the future contract
and currently fails at the schema-fixture gate.

## CLI Contract

`summary` requires two drug names and two exact SMILES strings. `batch --input`
must be a non-empty JSON list; each item must contain:

```text
drug1_name, drug2_name, drug1_smiles, drug2_smiles
```

Allowed optional fields are `cell_line`, `predicted_synergy`, `advanced`,
`disease_label`, and `model_fidelity_evidence`. Unknown fields and
`database`, `local_targets`, `drug1_targets`, or `drug2_targets` are rejected.
Missing output directories are created automatically.

After a reviewed fixture is production-approved, every formal run must also
provide independent raw PLATO export SHA-256 allowlist entries:

```powershell
python InterpreterPLATO/scripts/run_interpreter.py summary `
  --context-config InterpreterPLATO/data/production_context.json `
  --drug1 "Drug-A" `
  --drug2 "Drug-B" `
  --drug1-smiles "CCO" `
  --drug2-smiles "CCN" `
  --trusted-plato-sha256 "<drug-a-raw-export-sha256>" `
  --trusted-plato-sha256 "<drug-b-raw-export-sha256>"
```

Cache entries must use current cache schema v3 and pass an independently audited,
production-allowlisted r37 fixture gate, plus query-structure, release, raw-export
hash, reliable target, explicit-human, and protein-target gates. No fixture is
currently approved, so the present candidate cannot pass the first gate. Future
accepted outputs use
`result_schema_version=5.1` and
`credibility.formula_version=evidence_credibility_v4_verified_plato_only`.
Schema 3/4/5.0 and formula v1/v2/v3 files are historical audit artifacts only
and must not be compared directly with v5.1 scores.

## Checks

```powershell
python InterpreterPLATO/scripts/run_interpreter.py summary --help
python InterpreterPLATO/scripts/run_interpreter.py batch --help
python InterpreterPLATO/scripts/query_pair_and_explain.py --help
python InterpreterPLATO/scripts/query_plato_once.py --help
python InterpreterPLATO/scripts/import_plato_result.py --help
python -m pytest InterpreterPLATO/tests
```
