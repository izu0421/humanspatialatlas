#!/usr/bin/env Rscript
# Extract a Seurat .rds into plain files that Python can assemble into .h5ad.
#
# Deliberately uses NO Seurat/SeuratObject package. Those are not installed here, and more
# importantly an S4 method call on a Seurat v5 object requires the class definitions, which makes
# the conversion depend on matching package versions for every deposit. Reading the attributes
# directly needs only base R and Matrix, and works across v3, v4 and v5 layouts.
#
# Writes: counts.mtx (genes x cells), features.tsv, barcodes.tsv, obs.csv, coords.csv, info.json
suppressWarnings(suppressMessages(library(Matrix)))

args <- commandArgs(trailingOnly = TRUE)
f_in <- args[1]; d_out <- args[2]
dir.create(d_out, recursive = TRUE, showWarnings = FALSE)
note <- function(...) cat(..., "\n", sep = "")

x <- suppressWarnings(readRDS(f_in))
a <- attributes(x)
cls <- paste(class(x), collapse = ",")

# ---- counts ----------------------------------------------------------------
# v5 keeps matrices in assay$layers; v3/v4 in assay@counts or @data. Prefer raw counts.
# Only Seurat keeps assays as a plain named list. In a SummarizedExperiment the same slot is an
# S4 Assays object, and indexing it dispatches a method that needs the package installed -- which
# is exactly the dependency this script exists to avoid. Touch it only when it is a bare list.
pick_assay <- function(a) {
  tryCatch({
    if (is.null(a$assays) || isS4(a$assays) || !is.list(a$assays)) return(NULL)
    nm <- names(a$assays)
    if (is.null(nm)) return(NULL)
    pref <- c("Spatial", "RNA", "Nanostring", "Xenium", "SCT")
    hit <- nm[nm %in% pref]
    a$assays[[if (length(hit)) hit[1] else nm[1]]]
  }, error = function(e) NULL)
}
# Containers nest differently: Seurat keeps assays at the top level, SpatialExperiment hides the
# matrix under assays -> data (an S4 SimpleList) -> listData -> counts, and each S4 layer needs
# its own attributes() unwrap. Rather than hard-code a path per class, walk the object and take
# the first sparse matrix of plausible shape. Bounded depth so a self-referential slot cannot
# spin; the Seurat path above is still tried first because it knows which assay to prefer.
find_sparse <- function(obj, depth = 0, best = NULL) {
  if (depth > 6 || is.null(obj)) return(best)
  at <- attributes(obj)
  if (!is.null(at$Dim) && !is.null(at$x) && !is.null(at$i)) {      # dgCMatrix-like
    if (is.null(best) || prod(at$Dim) > prod(attributes(best)$Dim)) return(obj)
    return(best)
  }
  kids <- list()
  if (is.list(obj) && !isS4(obj)) kids <- obj
  if (!is.null(at$listData)) kids <- c(kids, at$listData)
  if (!is.null(at$data)) kids <- c(kids, list(at$data))
  if (!is.null(at$assays)) kids <- c(kids, list(at$assays))
  for (k in kids) best <- find_sparse(k, depth + 1, best)
  best
}

as1 <- pick_assay(a)
aa <- if (is.null(as1)) list() else attributes(as1)

M <- NULL; layer <- NA
if (!is.null(aa$layers)) {
  for (k in c("counts", "data", names(aa$layers))) {
    if (!is.null(aa$layers[[k]])) { M <- aa$layers[[k]]; layer <- k; break }
  }
} else {
  for (k in c("counts", "data")) {
    if (!is.null(aa[[k]]) && length(attributes(aa[[k]])$Dim)) { M <- aa[[k]]; layer <- k; break }
  }
}
if (is.null(M)) {                      # not a Seurat layout: search the container
  M <- find_sparse(x)
  layer <- "found"
}
if (is.null(M)) { note("ERROR no count matrix"); quit(status = 3) }
ma <- attributes(M)
dims <- ma$Dim
if (is.null(dims)) { note("ERROR matrix has no Dim"); quit(status = 4) }

# names: v5 holds them on the assay's cells/features LogMaps, v3/v4 on the matrix Dimnames
getnames <- function(obj) { d <- attributes(obj)$dimnames; if (is.null(d)) NULL else d[[1]] }
feats <- if (!is.null(ma$Dimnames) && length(ma$Dimnames[[1]])) ma$Dimnames[[1]] else getnames(aa$features)
cells <- if (!is.null(ma$Dimnames) && length(ma$Dimnames[[2]])) ma$Dimnames[[2]] else getnames(aa$cells)
# SummarizedExperiment keeps the names beside the matrix rather than on it
if (is.null(feats) && !is.null(a$elementMetadata)) {
  rn <- attributes(a$elementMetadata)$rownames
  if (is.character(rn) && length(rn) == dims[1]) feats <- rn
}
if (is.null(cells) && !is.null(a$colData)) {
  cn <- attributes(a$colData)$rownames
  if (is.character(cn) && length(cn) == dims[2]) cells <- cn
  if (is.null(cells)) {
    cl <- attributes(a$colData)$listData
    for (nm in c("_index", "barcode", "cell_id", "Barcode")) if (!is.null(cl[[nm]]) &&
        length(cl[[nm]]) == dims[2]) { cells <- as.character(cl[[nm]]); break }
  }
}
if (is.null(feats)) feats <- paste0("gene_", seq_len(dims[1]))   # positional, flagged below
if (is.null(cells)) cells <- paste0("cell_", seq_len(dims[2]))
# bin identifiers repeat in some SpatialExperiment deposits; R refuses non-unique row names and
# the Python side needs a usable index either way
cells <- make.unique(as.character(cells))
feats <- make.unique(as.character(feats))
if (is.null(feats) || is.null(cells)) { note("ERROR missing dimnames"); quit(status = 5) }

writeMM(M, file.path(d_out, "counts.mtx"))
writeLines(as.character(feats), file.path(d_out, "features.tsv"))
writeLines(as.character(cells), file.path(d_out, "barcodes.tsv"))

# ---- cell metadata ---------------------------------------------------------
md <- a$meta.data
if (is.null(md) && !is.null(a$colData)) {      # SummarizedExperiment-style
  cl <- attributes(a$colData)$listData
  keep <- cl[sapply(cl, function(z) is.atomic(z) && length(z) == dims[2])]
  if (length(keep)) {
    md <- as.data.frame(keep, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(md) <- cells
  }
}
if (!is.null(md)) {
  md <- as.data.frame(md)
  md <- data.frame(lapply(md, function(c_) if (is.list(c_)) sapply(c_, function(z) paste(unlist(z), collapse = ";")) else c_),
                   check.names = FALSE, stringsAsFactors = FALSE)
  write.csv(md, file.path(d_out, "obs.csv"), row.names = TRUE)
}

# ---- coordinates -----------------------------------------------------------
# Tried in order: v5 centroids, Visium image coordinates, then meta.data columns.
# A merged Seurat object carries ONE IMAGE PER SLIDE, each covering only its own cells, so the
# first image is not the sample's coordinates -- it is a fragment of them. Taking it alone gave
# 1,413 coordinates for a 3,580-cell object, which the Python side then rejected outright.
# Collect every image and stack them; the slide each came from is kept as a column so a merged
# object is recognisable downstream rather than silently flattened.
coords <- NULL; csrc <- "none"
if (!is.null(a$images) && length(a$images)) {
  parts <- list(); kinds <- character(0)
  for (nm in names(a$images)) {
    im <- attributes(a$images[[nm]]); got <- NULL
    if (!is.null(im$boundaries) && !is.null(im$boundaries$centroids)) {
      cc <- attributes(im$boundaries$centroids)
      if (!is.null(cc$coords)) {
        got <- as.data.frame(cc$coords)
        if (!is.null(cc$cells)) rownames(got) <- as.character(cc$cells)
        kinds <- c(kinds, "centroids")
      }
    }
    if (is.null(got) && !is.null(im$coordinates)) {        # Visium
      cd <- as.data.frame(im$coordinates)
      xy <- intersect(c("imagecol", "imagerow", "col", "row", "x", "y"), names(cd))
      if (length(xy) >= 2) {
        pick <- if (all(c("imagecol", "imagerow") %in% names(cd))) c("imagecol", "imagerow") else xy[1:2]
        got <- cd[, pick, drop = FALSE]
        kinds <- c(kinds, "coordinates")
      }
    }
    if (!is.null(got) && nrow(got)) {
      names(got)[1:2] <- c("x", "y")
      got$slide <- nm
      parts[[nm]] <- got[, c("x", "y", "slide")]
    }
  }
  if (length(parts)) {
    coords <- do.call(rbind, parts)
    rownames(coords) <- unlist(lapply(parts, rownames), use.names = FALSE)
    csrc <- paste0("images(", length(parts), "):", paste(unique(kinds), collapse = "/"))
  }
}
if (is.null(coords) && !is.null(a$int_colData)) {       # SpatialExperiment
  sc <- attributes(a$int_colData)$listData$spatialCoords
  if (!is.null(sc) && nrow(sc) == dims[2]) {
    coords <- as.data.frame(sc[, 1:2, drop = FALSE])
    names(coords) <- c("x", "y"); rownames(coords) <- cells
    csrc <- "int_colData:spatialCoords"
  }
}
if (is.null(coords) && !is.null(md)) {
  pats <- list(c("CenterX_global_px", "CenterY_global_px"), c("x_centroid", "y_centroid"),
               c("x_global_px", "y_global_px"), c("X", "Y"), c("x", "y"),
               c("imagecol", "imagerow"), c("sdimx", "sdimy"))
  for (p_ in pats) if (all(p_ %in% names(md))) {
    coords <- md[, p_, drop = FALSE]
    names(coords) <- c("x", "y")        # the Python side keys on lower-case x/y
    csrc <- paste0("meta.data:", paste(p_, collapse = ",")); break
  }
}
if (!is.null(coords)) {
  write.csv(coords, file.path(d_out, "coords.csv"), row.names = TRUE)
}

assay_nm <- tryCatch({
  n_ <- if (!isS4(a$assays) && is.list(a$assays)) names(a$assays)[1] else NA
  if (is.null(n_) || is.na(n_)) "unknown" else n_
}, error = function(e) "unknown")
cat(sprintf('{"class":"%s","assay":"%s","layer":"%s","n_genes":%d,"n_cells":%d,"coord_source":"%s","n_obs_meta":%d}\n',
            cls, assay_nm, layer, dims[1], dims[2], csrc,
            if (is.null(md)) 0 else nrow(md)),
    file = file.path(d_out, "info.json"))
note("OK ", dims[1], " genes x ", dims[2], " cells, coords=", csrc)
