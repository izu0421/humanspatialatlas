# HSA — to do

## Dashboard
- [x] **Headline line: "X cells across Y donors in Z tissues".** Rendered above the KPI tiles.
      Measured, human-only, loaded-samples-only: **126,833,893 single cells and 299,060,172
      spots/bins across ≥3,076 donors in 277 UBERON tissues** (5,633 samples, 1,111 datasets).
      Cells and spots are deliberately not summed: 259.8M of the units are Visium HD 2µm bins.
      The earlier 654M/2,637/289 figure in this file was wrong — it mixed species and counted
      samples whose matrix does not load. Donor count is a floor (61.7% of samples state one, and
      donor ids are study-local, so counted as (accession, donor_id) pairs).
- [x] Cells-per-technology panel ("Measured units per technology"), cells in blue and spots/bins
      in pink so the Visium HD bin contribution is visible rather than hidden in a total.
- [x] Load/QC status per technology ("What actually opens"), stacked as loads / no coordinates /
      no matrix found / read error. Reveals that Xenium's 37.9% load rate is a *coordinate*
      problem (716 of 978 failures), not a matrix one — a different fix from Visium's.

## Corpus
- [ ] Raise loader coverage. **Root cause found for the largest single block:** GEO re-gzips 10x
      parquet deposits, so Xenium centroids arrive as `cells.parquet.gz` and Visium HD positions as
      `tissue_positions.parquet.gz`. `_read_coord_table()` tested `p.suffix == ".parquet"`, which is
      `".gz"` for those, so a binary parquet went into `read_csv` and the exception was swallowed as
      `no_coords`. 827 such files are on disk (575 Xenium + 150 Visium HD + 102 others); 8 of 8
      previously-failing Xenium samples now load with 43k-253k cells each. Fixed; `qc_all(retry=
      "no_coords")` added to roll the fix out over exactly the affected samples.
      Still open afterwards: 3,039 `no_recognised_matrix` and whatever `no_coords` survives.
- [ ] Decide whether to materialise standardised `.h5ad` per sample (needs ~1–1.5 TB; /data at 97%).
- [ ] Re-run metadata harmonisation for datasets recovered after the archive fixes.
- [ ] Push refreshed dashboard + exports (one commit outstanding).

## Visium HD cell-level
- [ ] bin2cell: 13 of 19 eligible samples segmented (18/19 in progress, the 11 mm colon section).
      Fixed a silent skip: the output filename truncated at 120 characters, so three distinct
      lung-cancer post-Xenium samples collapsed onto one path and two were recorded as
      `already_done` against the first one's file. `_tag()` now appends an md5 of the full
      identity and migrates the existing files on next touch, so those two will be re-queued.
- [ ] 131 samples have an image but no 2 µm bins — would need `binned_outputs` re-download
      (8–12 GB each) to become eligible.
- [ ] 309 of 459 HD samples have no public full-res image at all; depositors ship the CytAssist
      image (~14–26 MB) and fiducial JPEGs, which cannot resolve nuclei.

## PPI benchmark
- [x] Atera benchmark running on two samples (breast, cervical) at a standardised 30,000-cell
      budget so every method is scored on the identical pair list.
- [x] LIANA+ bivariate added. It wins and it replicates: AUPRC 0.154 breast / 0.137 cervical
      against a 0.091 random floor. Specificity is its weakness (7.4% / 12.2% of impossible
      housekeeping x ligand pairs clear the decoy-95th threshold).
- [x] COMMOT added. It does not replicate: 0.113 breast but 0.083 cervical against a 0.091 floor,
      AUROC 0.453, i.e. slightly anti-correlated with CellPhoneDB in the second tissue.
- [ ] SpatialDM running. The earlier `IndexError` was ours: `globle_st_compute()` sizes its
      variance vector by counting pairs annotated `ECM-Receptor`/`Cell-Cell Contact`/
      `Secreted Signaling`, and we had labelled every pair `custom`, so `st` came back length-0.
- [ ] Celcomen still not evaluable at this scale: lr 1e-1 saturates every g2g entry to 1.0 and
      lr 1e-6 leaves it at initialisation. Report as "not evaluated", never as "Celcomen fails".
- [ ] Notebook with PR curves and the specificity plane (per-pair scores now persisted to
      `results/scores_{sample}.parquet`).
- [ ] Extend to a second platform once a multi-donor single-cell + wide-panel substrate exists.

## Known caveats to keep visible
- Public repo: earlier commits still contain `exports/hst_corpus_coverage.csv` with internal
  `/lustre` paths. History rewrite not done (needs explicit OK — irreversible).
- GeoMx (56 samples, 1,734 `.dcc.gz`) is region-level, not cell-level; currently flagged, not loaded.
