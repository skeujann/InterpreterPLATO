# Data And Provenance

`data/raw/` stores fixed-version source snapshots. `data/processed/context/`
stores the human PPI, expression, pathway, and disease context mapped to the
fixed HGNC release. These files support cancer-context reranking only. They do
not provide drug targets.

The active/default target source is one archived PLATO response per molecule.
The previous local demo target fixture has been moved to
`data/processed/archive/pre_plato_only/` for audit history.

```powershell
python InterpreterPLATO/scripts/build_context_artifacts.py all
python InterpreterPLATO/scripts/verify_manifest.py --scope runtime
```

Production context is loaded atomically from `data/production_context.json`.
Each resource should record source, version, `taxid=9606`, identifier namespace,
HGNC mapping version, file SHA-256, and content SHA-256. Runtime identity is
checked before explanation; mismatches fail closed. A source name alone is not
upstream authentication and must be checked against independent download records
or trusted manifests.

## PLATO Cache

Default location: `data/processed/plato_cache/`.

Current cache schema is `3`. It preserves immutable raw response bytes and
their SHA-256, plus normalized payload, normalized targets, query binding,
parser report, and audit metadata. Old cache schema v1/v2 entries are
archive-only and cannot be used for formal v5.1 output.

Archived files come from manual or explicitly user-triggered single-molecule
PLATO web queries. The interpreter itself does not fetch targets; the dedicated
helper permits one audited HTTPS form submission at a time. User attestation is
not provider authentication or schema verification.

A locally acquired Gemcitabine PLATO r37 target-fishing fixture may be placed
under `data/raw/plato/r37/golden_fixture/`. Raw responses and PDF reports are
intentionally excluded from this code-only repository. When supplied, the JSON
fixture verifies the machine-readable `chemblId / nome / score / validato /
simili` schema; it does not by itself prove every future query. Formal scoring
still requires query-structure binding and independent raw export SHA-256
allowlist entries. Do not handcraft JSON to bypass these gates or invent the
r37 schema.

## Scoped Locks And Provenance

- `runtime.lock.json`: default, fast verification scope. It covers the
  production config, its five runtime files, and the offline PLATO JSON fixture.
- `source.lock.json`: complete 35.45 GiB `data/raw` inventory. It is intentionally
  slow and is checked with `verify_manifest.py --scope source`.
- `manifest.lock.json`: detailed provenance and licensing companion retained for
  compatibility; it is not the default runtime verification scope.
- `manifest.example.json`: provenance template, not consumed at runtime.
- `manifest.local.json`: local migration-era snapshot; it cannot replace the
  locked manifest.

Templates describe context data only. They no longer include a local target
SQLite database. Rebuilding context artifacts changes QC timestamps and hashes;
review and update lock files intentionally, never just to satisfy a test.

## Result Boundary

Future formal results use `result_schema_version=5.1` and
`formula_version=evidence_credibility_v4_verified_plato_only`. Old schema
3/4/5.0, formula v1/v2/v3, cache v1/v2, and user-attested-only archives are not
valid v5.1 production inputs or score-comparison baselines.
