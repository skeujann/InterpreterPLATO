# Data bootstrap and offline rebuild

This document installs or reproduces the currently pinned snapshot. For an
upstream version upgrade, candidate build, biological QC review, lock update,
PLATO evidence refresh, and rollback, follow the Chinese
[data update operations manual](DATA_UPDATE_GUIDE.md).

InterpreterPLATO deliberately separates code, runtime context, and source data:

| Layer | Lock | Approximate size | Normal use |
|---|---|---:|---|
| Python code and tests | `requirements.lock` | small | install and CI |
| Runtime context plus offline parser-regression candidate fixture | `data/runtime.lock.json` | 7.64 MiB | default verification |
| Complete local `data/raw` snapshot | `data/source.lock.json` | 35.45 GiB | rebuild and provenance audit |
| ScoreRI cell-line registry | `data/scoreri_cell_lines.tsv` | 569 bytes | expression rebuild only |

The runtime bundle is not embedded in the Python wheel. The source bundle is
also not a Python dependency and should not be copied into a wheel or source
archive. Review `NOTICE` and every upstream license before distributing either
bundle.

The runtime lock proves byte identity, not fixture approval. The current PLATO
fixture manifest remains a candidate and the production schema allowlist is
empty, so formal scoring fails closed until an independently reviewed fixture is
approved in code, provenance, and tests together.

## 1. Install the code

Use CPython 3.9 through 3.12. From the `InterpreterPLATO` directory:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install --requirement requirements.lock
.venv\Scripts\python -m pip install --no-deps --editable .
. .\.venv\Scripts\Activate.ps1
python -c "import sys; print(sys.executable)"
```

If PowerShell policy blocks activation, keep the environment binding explicit by
replacing every later `python` with `.venv\Scripts\python` and every console
entry point with `.venv\Scripts\<entry-point>.exe`.

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
`--component context` or `--component target-mapping`, but the wrapper still
verifies all 12 registry inputs before either component. It writes directly to
the live `data/processed/context` directory and is therefore intended only for
reproducing the already pinned snapshot, not for upgrading releases.

Generated mapping gzip files use `mtime=0` and atomic replacement so identical
inputs produce identical bytes. QC reports intentionally contain audit times
and are not runtime trust anchors.

The rebuild never rewrites a lock automatically. If an intentional source or
builder change produces different output, stop using the wrapper and follow
`DATA_UPDATE_GUIDE.md`:

1. build into an inactive candidate directory;
2. review build QC, rejection sidecars, biological coverage, source, and code;
3. promote a versioned runtime bundle;
4. update production configuration and all affected locks in the documented order;
5. run runtime/source verification and the full test suite.

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
