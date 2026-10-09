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
- [x] Two Atera tumours (breast, cervical) at a standardised 30,000-cell budget, all methods on an
      identical pair list.
- [x] **Answer: the two leading methods are the same statistic.** SpatialDM 0.1550 / 0.1483 and
      LIANA+ bivariate 0.1540 / 0.1369 against a 0.091 floor — indistinguishable on breast
      (Δ 0.001). Both are bivariate Moran's R; the kernel choice does not matter. Everything
      structurally different is at or below the floor: COMMOT flips (0.1125 → 0.0826, AUROC 0.453),
      naive and residualised cross-covariance sit at the floor. Ceiling of the field ≈ 1.6-1.7x
      random while calling 6-13% of pairs that cannot physically interact.
- [x] SpatialDM fixed. The `IndexError` was ours: `globle_st_compute()` sizes its variance vector by
      counting pairs annotated `ECM-Receptor`/`Cell-Cell Contact`/`Secreted Signaling`, and we had
      labelled every pair `custom`, so `st` came back length-0. The label also selects the kernel.
- [x] **Negative control was abundance-confounded — corrected.** The 20 housekeeping genes have mean
      detection 0.152 vs 0.066 for positives. Naive cross-covariance's "22% of impossible pairs
      called", which we read as the co-location confound and as justification for residualising,
      drops to **5.7%** against closed-compartment partners drawn from the same detection decile;
      co-location's rate *doubles* (3.1 → 6.0%). So the naive statistic is abundance-biased, not
      architecture-dominated, and residualisation bought little real specificity while costing the
      discriminative signal. New class: 1,740 GO closed-compartment genes minus anything ever
      annotated surface-exposed or secreted (`ppib/truth.py: intracellular_genes()`).
- [x] Seed variance measured (`run_noise_floor.py`): AUPRC moves 0.012-0.019 between seeds for the
      clustering-dependent methods, 0.001 for naive cross-cov. **Nothing under ~0.02 AUPRC is a
      result.** Two earlier "independent replicates" shared `seed=0` and so were near-duplicates —
      fixing the seed for cross-method comparability is not reproducibility.
- [ ] Corrected (prevalence-matched) FDR for SpatialDM, LIANA+ and COMMOT — running. The 7.4-12.7%
      figures currently quoted for them are the unmatched ones and will move.
- [ ] Replicated depth sweep (4 budgets x 3 seeds x 2 tumours). Needed because residualised
      cross-cov gave 0.112 at 80k and 0.097 at 30k, but its seed spread alone covers most of that.
- [ ] Celcomen still not evaluated: saturates at the tutorial lr, sits at initialisation below it,
      and costs ~5 h (44 s/epoch x 200 x 2 samples). Report as "not evaluated", never "fails".
- [ ] Execute `notebooks/ppi_benchmark.ipynb` (20 cells, 8 figures) once the above land.
- [ ] Possible upgrade: build the negative class from Zhang et al. 2026 (aem7299) validated
      localisation markers — 2,011 organelle + 3,587 topology — instead of GO terms.

## Known caveats to keep visible
- Public repo: earlier commits still contain `exports/hst_corpus_coverage.csv` with internal
  `/lustre` paths. History rewrite not done (needs explicit OK — irreversible).
- GeoMx (56 samples, 1,734 `.dcc.gz`) is region-level, not cell-level; currently flagged, not loaded.
