<img src="logo.svg" width="96" align="right" alt="HSA logo">

# HSA — Human Spatial Atlas

A uniformly collected catalogue of public **human** spatially-resolved transcriptomics data
(Visium / Visium HD, Xenium, CosMx, MERFISH, Stereo-seq, Slide-seq, …), modelled on
[scBaseCount](https://doi.org/10.1101/2025.02.27.640494). An LLM agent (Claude Sonnet 5.5) finds and
curates datasets; everything that touches disk is deterministic code.

```
 sources ──► candidates ──► curator agent ──► deterministic gate ──► catalogue ──► downloader ──► data/
              (queue)        (1 per dataset)   (URL was listed?        (SQLite)    (matrix + coords
 GEO      : E-utilities keyword union                no images/transcripts?                only; remote-zip
 lists    : reference dataset lists                  'ready' needs both?)                  member reads;
 scouts   : 10x, Zenodo, HuBMAP, HTAN, vendors…                                             tar streaming)
 CELLxGENE: curation API (no LLM) ─────────────────────────────────────────►
```

Current phase: **cell/spot-by-gene count matrix + spatial coordinates per sample**. Images, raw reads
and per-transcript tables are linked in the catalogue but not downloaded.

## Outputs

| file | content |
|---|---|
| `exports/hsa_datasets.csv` | one row per dataset: source, accession, technology, tissue, disease, verdict, samples, files, GB, links |
| `exports/hsa_files.csv` | one row per file: sample, URL, role (matrix / coords / matrix+coords / bundle), download status |
| `exports/hsa_sources.csv` | database registry: URL, access method, estimated # human datasets |
| `exports/hsa_samples.csv` | one row per sample: harmonised species, tissue (UBERON), sex, age / development stage (HsapDv), disease levels 1-3 (MONDO), ethnicity (HANCESTRO), donor, provenance |
| `index.html` | self-contained dashboard, served by GitHub Pages — https://www.yizhouyu.com/humanspatialatlas/ |

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
