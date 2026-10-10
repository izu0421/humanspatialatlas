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
- [x] **Gzipped-parquet coordinate fix rolled out.** GEO re-gzips 10x parquet deposits, so Xenium
      centroids arrive as `cells.parquet.gz` and Visium HD positions as `tissue_positions.parquet.gz`;
      `_read_coord_table()` tested `p.suffix == ".parquet"`, which is `".gz"`, so a binary parquet
      went to `read_csv` and the exception was swallowed as `no_coords`. Retried all 1,894 affected
      samples: **+733 now load (57.0% -> 63.3%)**, Xenium's loaded count doubled 596 -> 1,216
      (+72,638,400 single cells, +57% on the atlas total), Visium HD +92 samples (+151M bins).
- [x] Three more matrix formats read (2026-10-10): gzipped 10x `.h5` (h5py cannot read a gzip
      stream, so it is decompressed first), Xenium `cell_feature_matrix.zarr.zip` (a zipped zarr
      store holding the matrix as CSR-over-features and/or CSC-over-cells; take whichever is
      present and transpose, obs_names positional so coordinates join by order), and mtx trios
      named `counts.mtx` / `count_matrix_sparse.mtx` with an orientation chosen from the sidecar
      lengths rather than assumed. **+45 samples, +8.5M cells.**
      Estimating from file-extension counts predicted ~280 samples and delivered 45 - extension
      counts are not sample counts, and several of those samples had other problems too. Count
      samples matching the new pattern before writing the reader next time.
- [ ] Remaining: 2,899 `no_recognised_matrix`, of which 1,640 are GeoMx (region-level by design,
      out of scope) and ~346 are R `.rds` Seurat objects that would need rpy2. ~900 unclassified.
      Plus ~1,040 residual `no_coords` in Visium `tissue_positions` variants and CosMx conventions.
      A deposit whose features.tsv does not describe its matrix (GSE287459: 37,082-row matrix,
      18,085-line features) now fails with that stated, rather than an opaque length error.
- [ ] Decide whether to materialise standardised `.h5ad` per sample (needs ~1–1.5 TB; /data at 97%).
- [ ] Re-run metadata harmonisation for datasets recovered after the archive fixes.
- [ ] Push refreshed dashboard + exports (one commit outstanding).

## Visium HD / Xenium uniform preprocessing
- [x] **Panel registry** (`hsa/panels.py`). Gene names read without the matrix (feature_type-filtered
      h5, features.tsv, h5ad var index): 1,536 samples in 22 s vs ~7 h via `load_sample`, validated
      exact against the loader on 5 of 6 samples (6th differs only by its duplicate renaming).
      **1,536 samples -> 164 distinct gene panels.** Panel SIZE is not panel identity: the 164
      samples reporting 541 genes span four panels sharing as little as 11% of their genes. Also
      found `n_genes` in sample_qc.csv is inconsistent for 446 samples because Xenium/CellRanger h5
      files carry control codewords alongside real genes (one 541-feature file: 300 Gene Expression
      + 241 controls); the registry stores control-stripped counts.
- [x] **No pan-Xenium gene space exists.** Intersection of all 34 panels with >=300 genes is ONE
      gene; reaching ~300 shared genes requires discarding more than half the panels. Bears directly
      on TERRA/Nicheformer-style shared-gene-space training.
- [x] **Xenium cohorts**: 5K family 248 samples / 46 datasets / 4,277 genes; 377-family 231 / 43 /
      377 genes; plus 12 smaller cohorts. See `runs/panel_cohorts.csv`.
- [x] **Visium HD bin size resolved exactly.** Only 96 of 320 samples declare one. Geometry cannot
      settle it (per-sample px/um ratio 0.25-6.17; three pitch estimators all plateaued near 57%),
      but the grid is deterministic: Space Ranger lays a fixed bin count over the capture area, so
      array_row's extent is an exact key -- 3350/838/419 rows for 2/8/16 um on the 6.5 mm slide,
      5500/1375/688 on the 11 mm slide. 100% agreement on all 94 declared samples.
- [x] **Only 205 of 320 samples labelled "Visium HD" are Visium HD** (`reclassify()`): 64 are on
      classic Visium's 78x128 spot grid, 33 are Xenium (continuous coordinates from a
      cells.parquet), 11 are cell-level WTA, 6 cell-level unknown, 1 unresolved. The dashboard's
      HD counts need to use the resolved technology, not the depositor's label.
- [x] **Uniform HD cohorts** (`runs/hd_cohorts.csv`): 2 um 69 samples / 40 datasets / 509.1M bins;
      8 um 96 / 43 / 33.3M; 16 um 40 / 30 / 5.2M. All three share the same 17,699 genes, and all
      205 genuine HD samples are covered with none excluded.
- [ ] Wire the resolved technology into `export.py` so the dashboard stops counting 64 Visium and
      33 Xenium samples as Visium HD (~21M units mis-attributed; the cell/bin split also shifts).
- [x] Coordinate scale resolved per sample (`panels.coordinate_scale()`, runs/coordinate_scale.csv).
      um_per_unit = (grid_rows x bin_um) / coordinate extent. Among the 205 genuine HD samples,
      187 are in full-resolution pixels (median 0.287 um/unit) and 18 already in microns (median
      0.95, which is the check that the derivation is sound rather than rescaling noise). The
      earlier "225 pixels / 92 microns" figure was computed over all 320 labelled samples with an
      extent heuristic and is superseded.
- [ ] Apply the scale when materialising, rather than mutating stored files.
- [ ] Materialise the cohorts as standardised .h5ad. Blocked on disk: /data at 98%, 1.4 TB free,
      full corpus estimated 1-1.5 TB -- plan is the three HD and two large Xenium cohorts only,
      int32 CSR, no dense layers.
- [ ] 2 samples on an unrecognised 280-row grid (GSE325706) left flagged, not guessed.

## Visium HD cell-level
- [x] GEX rasterisation resolution: `gex_mpp` default changed from bin2cell's 2.0 to **0.5**.
      At 2 um/px a 10 um nucleus spans 5 px and StarDist merges them. Measured on GSM9937429
      (9.83M bins): mpp 2.0/1.0/0.5 -> 9,829 / 17,959 / 40,267 nuclei, against 74,928 cells for
      the SAME sample segmented from H&E once its image finally downloaded. So 2.0 recovers 13%
      of the H&E yield and 0.5 about 54%. Costs 16x the pixels (~1 min -> ~10 min per sample).
      Note cells-per-bin is NOT a platform constant and cannot be used to detect over-segmentation
      across tissues: H&E rates range 0.70% (GSE342738 brain) to 4.93% (10x breast cancer).

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
