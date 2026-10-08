# HSA — to do

## Dashboard
- [ ] **Headline line: "X cells across Y donors in Z tissues".** Put it at the top of the page,
      above the KPI tiles, as the one sentence that says what HSA is. Current values:
      **654,095,070 cells · 2,637 donors · 289 UBERON tissues** (also available: 6,643 samples with
      cell counts, 379 MONDO diseases, 2,102 datasets).
      Needs `runs/sample_qc.csv` joined into `hsa/export.py`; donor count comes from
      `sample_meta.donor_id` (only ~24% of samples state one, so the honest phrasing is
      "≥2,637 donors" or report samples alongside).
- [ ] Add a cells-per-technology panel (Visium HD bins vs true single cells — label the difference,
      since 418M of the 654M are 2/8/16 µm bins, not cells).
- [ ] Show load/QC status per technology so a user can see what actually opens.

## Corpus
- [ ] Raise loader coverage: 57% of samples load cleanly. Remaining failures are
      3,037 `no_recognised_matrix` and 1,734 `no_coords` — more format families, not bad data.
- [ ] Decide whether to materialise standardised `.h5ad` per sample (needs ~1–1.5 TB; /data at 97%).
- [ ] Re-run metadata harmonisation for datasets recovered after the archive fixes.
- [ ] Push refreshed dashboard + exports (one commit outstanding).

## Visium HD cell-level
- [ ] bin2cell on the 19 samples that have image + 2 µm bins (image reading just fixed:
      gzipped TIFF and BigTIFF both return empty from `cv2.imread` rather than raising).
- [ ] 131 samples have an image but no 2 µm bins — would need `binned_outputs` re-download
      (8–12 GB each) to become eligible.
- [ ] 309 of 459 HD samples have no public full-res image at all; depositors ship the CytAssist
      image (~14–26 MB) and fiducial JPEGs, which cannot resolve nuclei.

## PPI benchmark
- [ ] Finish the Atera benchmark (co-location baseline AUPRC 0.091 vs 0.048 random; cross-covariance
      methods still scoring).
- [ ] Add LIANA+ bivariate and SpatialDM as published comparators.
- [ ] Extend to a second platform once a multi-donor single-cell + wide-panel substrate exists.

## Known caveats to keep visible
- Public repo: earlier commits still contain `exports/hst_corpus_coverage.csv` with internal
  `/lustre` paths. History rewrite not done (needs explicit OK — irreversible).
- GeoMx (56 samples, 1,734 `.dcc.gz`) is region-level, not cell-level; currently flagged, not loaded.
