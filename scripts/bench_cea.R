#!/usr/bin/env Rscript

# Benchmark terra-based CEA vs matrix-based CEA across grid sizes
# Run from project root: Rscript scripts/bench_cea.R

suppressPackageStartupMessages({
  library(terra)
})

# Load package code without installation (R scripts only)
r_scripts <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)
invisible(lapply(r_scripts, source, chdir = TRUE))

set.seed(123)

# Matrix-based CEA (in-memory)
cea_matrix <- function(drivers, vc, sensitivity) {
  dr_mat <- values(drivers, mat = TRUE) # ncell x ndr
  vc_mat <- values(vc, mat = TRUE) # ncell x nvc
  ce_total <- numeric(nrow(dr_mat))
  for (j in seq_len(ncol(dr_mat))) {
    exp_j <- vc_mat * dr_mat[, j] # exposure for driver j (element-wise per vc)
    eff_j <- sweep(exp_j, MARGIN = 2, sensitivity[, j], `*`) # apply sensitivity
    ce_total <- ce_total + rowSums(eff_j, na.rm = TRUE)
  }
  ce_total
}

# Terra-based CEA (using package functions)
cea_terra <- function(drivers, vc, sensitivity, cores = NULL, filename = NULL, engine = "terra") {
  ce <- cea(drivers, vc, sensitivity, exportAs = "SpatRaster", cores = cores, filename = filename, engine = engine)
  cea_extract(ce, cumul_fun = "full")
}

run_case <- function(side, ndr, nvc) {
  message(sprintf("Grid %dx%d, drivers=%d, vc=%d", side, side, ndr, nvc))
  tmpl <- rast(ncols = side, nrows = side, xmin = 0, xmax = side, ymin = 0, ymax = side, crs = "EPSG:3857")
  drs <- do.call(c, lapply(seq_len(ndr), function(i) { r <- tmpl; values(r) <- runif(ncell(r)); names(r) <- paste0("driver", i); r }))
  vcs <- do.call(c, lapply(seq_len(nvc), function(i) { r <- tmpl; values(r) <- runif(ncell(r)); names(r) <- paste0("vc", i); r }))
  sens <- matrix(runif(nvc * ndr), nrow = nvc, ncol = ndr, dimnames = list(names(vcs), names(drs)))

  tm_mat <- try(system.time(res_mat <- cea_matrix(drs, vcs, sens))["elapsed"], silent = TRUE)
  tm_terra <- try(system.time(res_terra <- cea_terra(drs, vcs, sens, engine = "terra"))["elapsed"], silent = TRUE)

  # Only run parallel/file-backed if terra succeeded
  tm_terra_par <- NA
  if (!inherits(tm_terra, "try-error")) {
    cores_use <- min(4, max(1, parallel::detectCores() - 1))
    tmpfile <- tempfile(fileext = ".tif")
    terraOptions(parallel = TRUE, memfrac = 0.8)
    tm_terra_par <- try(system.time(
      res_terra_par <- cea_terra(drs, vcs, sens, cores = cores_use, filename = tmpfile, engine = "terra")
    )["elapsed"], silent = TRUE)
    unlink(tmpfile)
  }

  list(
    tm_mat = tm_mat,
    tm_terra = tm_terra,
    tm_terra_par = tm_terra_par
  )
}

# Sweep sizes to see where matrix fails or slows
grid_tests <- list(
  list(side = 500, ndr = 5, nvc = 5),
  list(side = 1000, ndr = 10, nvc = 10),
  list(side = 2000, ndr = 10, nvc = 10)
)

results <- lapply(grid_tests, function(cfg) do.call(run_case, cfg))

for (i in seq_along(grid_tests)) {
  cfg <- grid_tests[[i]]
  res <- results[[i]]
  msg <- sprintf("Grid %dx%d d=%d vc=%d | matrix: %s | terra: %s | terra+par: %s",
                 cfg$side, cfg$side, cfg$ndr, cfg$nvc,
                 if (inherits(res$tm_mat, "try-error")) "FAILED" else sprintf("%.3f", res$tm_mat),
                 if (inherits(res$tm_terra, "try-error")) paste("FAILED:", conditionMessage(attr(res$tm_terra, "condition"))) else sprintf("%.3f", res$tm_terra),
                 if (inherits(res$tm_terra_par, "try-error") || is.na(res$tm_terra_par)) "FAILED" else sprintf("%.3f", res$tm_terra_par))
  message(msg)
}
