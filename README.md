<img src="logo.svg" width="96" align="right" alt="HSA logo">

# HSA — the Human Spatial Atlas and the Human Spatial Agents

Two things, built together.

**The Human Spatial Atlas** is a uniformly processed catalogue of public **human** spatially-resolved
transcriptomics (Visium / Visium HD, Xenium, CosMx, MERFISH, Stereo-seq, Slide-seq, GeoMx, …):
cell- or spot-by-gene matrices with spatial coordinates and harmonised metadata, one record per sample.

**The Human Spatial Agents** are the pipeline that builds it and keeps it current. They scout the
repositories for new deposits, judge what is usable, fetch it, read it into one shape, harmonise the
labels against public ontologies, and derive what the depositor did not state. Every stage is
re-runnable, so the Atlas tracks the literature rather than being a snapshot of one afternoon.

Modelled on [scBaseCount](https://doi.org/10.1101/2025.02.27.640494). The agent (Claude Sonnet 5.5)
decides *what* to take; everything that touches disk is deterministic code.

```
          ┌──────────────────── the Human Spatial Agents ────────────────────┐
 SCOUT    13 lanes: GEO, Zenodo, figshare, Mendeley, CELLxGENE, 10x, HuBMAP,
          HTAN, ArrayExpress, HCA, Broad SCP, vendor showcases, aggregators
             │
 CURATE   one judgement per candidate: human? spatial? matrix + coordinates
          openly downloadable?  → deterministic gate → catalogue (SQLite)
             │
 FETCH    matrix and coordinate files per sample; resume, throughput floor,
          remote-zip member reads, tar streaming, archives unpacked
             │
 STANDARDISE  ten matrix formats → one shape; coordinates found however the
          platform stores them; every sample QC'd so failures are visible
             │
 HARMONISE  tissue · disease · sex · age · development stage · ethnicity
          → UBERON, MONDO, HsapDv, HANCESTRO (OLS4); sex also inferred from
          expression where the gene panel allows it
             │
 DERIVE    gene-panel identity · bin size from the bin grid · coordinate
          scale in microns · platform corrected from file evidence ·
          Visium HD segmented from 2 µm bins into cells (bin2cell)
          └──────────────────────────────┬──────────────────────────────────┘
                                         ▼
                            the Human Spatial Atlas
```

Current phase: **cell/spot-by-gene count matrix + spatial coordinates per sample**. Images, raw reads
and per-transcript tables are linked in the catalogue but not downloaded.

## Keeping it current

The agents are re-runnable, and one cycle is scheduled weekly:

```
17 3 * * 1   runs/hsa_cycle.sh      # python run_hsa.py cycle --budget 40
```

A cycle is scout → curate → download → QC → panels → materialise → export, each stage
resumable and writing its own ledger, logged to `runs/cycles/cycle_<stamp>.log`. Two guards
make it safe to leave unattended: the agent stages stop at both a per-cycle budget and an
absolute lifetime cap, so a looping lane cannot drain the account; and materialisation
refuses to start below a disk floor, because deleting raw files to make room for their
replacements is only safe while there is room for the replacements.

**What it does not do.** It acquires and judges autonomously, but it does not diagnose or
repair. Every matrix format it reads was added by hand, and a failing sample is recorded
faithfully and then stays failed. Worse, it cannot notice a result that is wrong but
plausible: segmenting Visium HD at 2 µm per pixel produced valid files with eight times
too few cells, and that surfaced only by comparing against H&E siblings. Treat the load
rate and the derived counts as things to check, not things the system defends.

## Outputs

Catalogue tables are generated locally by `python -m hsa.export` and are not published in this repository.

| file | content |
|---|---|
| `exports/hsa_master_sheet.csv` | **the master sheet** — one row per sample, 74 columns: identity, platform (as labelled and as resolved from file evidence), processing state (`usable` / `standardised` / `cell_called`), scale, harmonised biology with ontology ids, derived technical facts (panel id, bin size, capture area, µm per coordinate unit, inferred sex), artefact paths, and study context. Rebuilt on every `python -m hsa.export`. |
| `exports/hsa_datasets.csv` | one row per dataset: source, accession, technology, tissue, disease, verdict, samples, files, GB, links |
| `exports/hsa_files.csv` | one row per file: sample, URL, role (matrix / coords / matrix+coords / bundle), download status |
| `exports/hsa_sources.csv` | database registry: URL, access method, estimated # human datasets |
| `exports/hsa_samples.csv` | one row per sample: harmonised species, tissue (UBERON), sex, age / development stage (HsapDv), disease levels 1-3 (MONDO), ethnicity (HANCESTRO), donor, provenance |
| `index.html` | password-protected dashboard, served by GitHub Pages — https://www.yizhouyu.com/humanspatialatlas/ |

Verdicts: `ready` · `bundle_only` (matrix/coords only inside an archive; extracted by `bundles`) ·
`restricted` (login / token / registration needed) · `no_processed_data` · `not_spatial` · `not_human`.

## Run

```bash
micromamba env create -f environment.yml
echo "sk-ant-..." > ../api.txt                     # Anthropic key, read by hsa/config.py (never commit it)
python run_hsa.py harvest                          # CELLxGENE (deterministic)
python -c "from hsa import harvest_geo; harvest_geo.run()"    # queue GEO candidates
python -c "from hsa import harvest_xlsx; harvest_xlsx.run()"  # queue a reference dataset list (xlsx)
python run_hsa.py curate --workers 12 --budget 50  # Sonnet curator, global $ cap
python run_hsa.py scout                            # Sonnet scouts for the other databases
python run_hsa.py download                         # matrix + coords files
python run_hsa.py bundles                          # extract matrix + coords members from tar/zip archives
python run_hsa.py status
python -m hsa.metadata run                         # per-sample metadata harmonisation
python -m hsa.export                               # CSVs + dashboard
```

## Design notes

- **LLM decides, code acts.** The curator only records URLs that a listing tool actually returned; the
  gate rejects images, transcript tables and segmentation files, and a `ready` verdict needs both a matrix
  and coordinates. Rejections go back to the agent as tool errors.
- **Remote zips are never downloaded whole.** Xenium/CosMx `outs.zip` bundles (often >10 GB) are read
  via HTTP range requests; only `cell_feature_matrix.h5` / `cells.*` members are extracted.
- **Blocked sites.** 10xgenomics.com sits behind a bot checkpoint; the curator reads the Internet Archive
  snapshot of the dataset page to get exact file URLs, then HEAD-checks them live.
- **Not reachable without credentials:** Vizgen MERSCOPE showcase (Google sign-in), Dryad (API token),
  Google Drive folders, HTAN controlled-access levels.
- Cost so far: ≈ $0.035 per GEO series curated with Sonnet 5.5.
