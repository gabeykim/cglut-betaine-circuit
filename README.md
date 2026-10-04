# Candidate betaine responsive circuits in C. glutamicum

This computational study supplies 97 candidate DNA fragments and the analysis
that led to their selection: 73 output fragments and 24 sensor candidates.
The 95-record baseline became 97 through two added operator-scramble controls.
No functioning circuit, engineered strain or experimental induction is reported.

This is the public computational resource, version 0.6.0, at
https://github.com/gabeykim/cglut-betaine-circuit. Archived release: [version 0.6.0 on Zenodo](https://zenodo.org/records/23132170), DOI **10.5281/zenodo.23132170**.
The manuscript remains a draft. See
DEPOSIT_STEPS.txt inside the complete source archive for the remaining publication steps and declaration checks.

## Download the complete resource

The complete 427-file resource is in `cglut-betaine-circuit-v0.6.0.tar.xz`.
Download and extract it before following the full reproduction instructions:

```bash
tar -xJf cglut-betaine-circuit-v0.6.0.tar.xz
cd cglut-betaine-circuit
```

Selected current candidates and audit scripts are also shown directly in this
repository for browsing. Historical project inputs, outputs, logs, figures and
the draft manuscript are retained in the archive. Checksums in the extracted
`release_sha256.csv` cover the archive contents.

## What is included

- `candidates/`: CSV/FASTA, component maps, baseline, edits and design rationales.
- `project/`: original analysis scripts, retained source inputs and historical outputs.
- `publication_code/`: portable candidate audit, rerun driver, P1 builder and figures.
- `reproduction/`: environment pins, logs, SHA-256 comparisons and label-correction patch.
- `independent_audit/`: separately implemented checks of sequence/component consistency.
- `figures/`: figure code outputs and underlying numerical values.
- `manuscript_draft.txt`: current text; journal declarations remain unfinished.

Historical reports retain earlier interpretations. Read
`archival_interpretation_corrections.txt` before using them as scientific claims.

## Run the candidate audit

Using Python 3.12, from this directory:

```bash
python publication_code/audit_candidates.py candidates ../cglut_candidate_audit
```

This uses only the standard library. It checks sequence identities, lengths,
components and documented edits. It does not test biological function.

## Reproduce the selected computational stages

The recorded environment was Linux, Python 3.12.14 and HMMER 3.4. Install the
exact Python package versions in `requirements.txt` in an isolated environment:

```bash
python3.12 -m venv ../cglut_env
source ../cglut_env/bin/activate
python -m pip install -r requirements.txt
hmmscan -h
```

Some cached third-party inputs are deliberately excluded from this public
resource. Read `THIRD_PARTY_NOTICES.txt` and `excluded_inputs.csv`.
If you hold the complete review archive lawfully, restore the exact inputs:

```bash
python publication_code/restore_review_inputs.py /path/to/Glutamicum_Preprint_Draft
python publication_code/reproduce.py /absolute/path/to/a/new/run
```

The destination must not already exist. Current source downloads can differ
from historical inputs and are not interchangeable for byte-level comparisons.
The complete cached-input package produced 38/38 selected matching files,
including 5/5 checklist files. This was AI-executed; personal author verification
is pending. No source refetch, DeepTMHMM inference or edited-P1 folding is implied.

## Figures

Install `requirements-figures.txt`, then run
`python publication_code/make_figures.py`. The scientific target is the plotted
numeric data; binary image metadata can vary across rendering environments.

## Reuse

Original code: MIT. Original design/analysis data: CC BY 4.0, subject to source
rights. Consult LICENSE-SCOPE.txt and THIRD_PARTY_NOTICES.txt. The manuscript
publication license is still pending. Contact: gabe.y.kim@gmail.com.
