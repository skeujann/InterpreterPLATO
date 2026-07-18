# Data bootstrap and offline rebuild

InterpreterPLATO deliberately separates code, runtime context, and source data:

| Layer | Lock | Approximate size | Normal use |
|---|---|---:|---|
| Python code and tests | `requirements.lock` | small | install and CI |
| Production runtime context plus JSON fixture | `data/runtime.lock.json` | 7.64 MiB | default verification |
| Complete local `data/raw` snapshot | `data/source.lock.json` | 35.45 GiB | rebuild and provenance audit |
| ScoreRI cell-line registry | `data/scoreri_cell_lines.tsv` | 569 bytes | expression rebuild only |

The runtime bundle is not embedded in the Python wheel. The source bundle is
also not a Python dependency and should not be copied into a wheel or source
archive. Review `NOTICE` and every upstream license before distributing either
bundle.

## 1. Install the code

Use CPython 3.9 through 3.12. From the `InterpreterPLATO` directory:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install --requirement requirements.lock
.venv\Scripts\python -m pip install --no-deps --editable .
```

`requirements.lock` is an exact version lock recorded for Windows x86-64. It
does not claim cross-platform wheel hashes. Re-resolve and review the lock when
the OS, architecture, or Python minor version changes.

## 2. Plan acquisition

The registry is small, reviewable, and contains expected paths, sizes, SHA-256
digests, official landing pages, and redistribution warnings:

```powershell
plato-bootstrap-data plan
```

The command does not access the network. Acquire each file from its listed
official release page, respecting access controls and license terms, then place
it at the expected path. Do not store credentials, signed download URLs, or
session cookies in the registry or a lock file.

COSMIC is license-gated and must be obtained directly by an authorized user.
ChEMBL attribution files must accompany any permitted redistribution. DepMap,
STRING, Reactome, HGNC, MSigDB, GDSC, PRISM, OncoKB, GO, and PLATO artifacts
retain their own upstream terms; this project does not relicense them.

Some official downloads are containers rather than the exact registered input:
extract the ChEMBL SQLite archive to the registered `chembl_37.db` path, extract
the Reactome GMT member from its versioned ZIP, and decompress the authorized
COSMIC payload. Verification always applies to the final registered bytes.

## 3. ScoreRI cell-line registry

The expression build needs only the 81 requested ScoreRI cell-line names. They
are pinned in the package-local identifier registry:

```text
data/scoreri_cell_lines.tsv
```

The registry was derived from the sorted unique `cell_line_name` values in the
historical `ScoreRI/data/combo_data.tsv`; its source digest is recorded in
`data/source.lock.json`. The 26 MB combination table is no longer a build
dependency, so an independent checkout needs no sibling ScoreRI repository.
Verify the small registry with the source scope:

```powershell
python scripts\verify_manifest.py --scope source
```

## 4. Verify locks

Fast runtime verification is the default:

```powershell
plato-verify-data
# equivalent:
python scripts\verify_manifest.py --scope runtime
```

The runtime scope covers `data/production_context.json`, PPI, expression,
Reactome pathways, the disease set, the ChEMBL-to-HGNC target mapping, and the
two small PLATO fixture files required by offline tests.

The source check hashes the entire local raw snapshot and can take substantial
time and I/O:

```powershell
python scripts\verify_manifest.py --scope source
plato-bootstrap-data verify
```

The bootstrap verifier checks the rebuild-required subset from
`data/source_registry.json`; the source lock additionally covers every retained
file in `data/raw`.

## 5. Rebuild runtime context offline

After every registered input passes size and SHA-256 checks:

```powershell
plato-bootstrap-data rebuild `
  --component all
```

The wrapper calls the existing context and target-mapping builders with every
input path spelled out. It performs no download and finishes by running the
runtime lock verifier. Components can be rebuilt separately with
`--component context` or `--component target-mapping`.

Generated mapping gzip files use `mtime=0` and atomic replacement so identical
inputs produce identical bytes. QC reports intentionally contain audit times
and are not runtime trust anchors.

The rebuild never rewrites a lock automatically. If an intentional source or
builder change produces different output:

1. review build QC and biological coverage changes;
2. review the exact source and code revision;
3. update `data/production_context.json` and `data/runtime.lock.json` together;
4. run `python scripts\verify_manifest.py --scope runtime`;
5. commit code, configs, locks, and the experiment record as one reviewed unit.

Do not edit hashes merely to silence a failed check.

## 6. SMILES inputs are not bundled

No drug-name-to-SMILES table ships with InterpreterPLATO, and the pair helper
does not scan ambient sibling projects. Every run must either supply both
SMILES directly or pass one or more explicit tables:

```powershell
plato-query-pair `
  --drug1 Drug-A --drug2 Drug-B `
  --smiles-lookup C:\datasets\drug_smiles.csv `
  --resolve-only
```
