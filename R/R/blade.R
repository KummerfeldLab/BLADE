#' Read standard Visium inputs
#'
#' @param dir Path to a Spaceranger output folder containing filtered_feature_bc_matrix/ and spatial/.
#' @param h5_path Optional path to filtered_feature_bc_matrix.h5.
#' @param neighbor_mode "radius" or "k" for neighbor selection.
#' @param k_neighbors Number of neighbors for k mode.
#' @param neighbor_radius Optional radius for neighbor search.
#' @return A list with matrix, positions, inclusion mask, and neighbor index.
#' @export
blade_read_visium <- function(dir, h5_path = NULL, neighbor_mode = "radius", k_neighbors = 6, neighbor_radius = NULL) {
  matrix_path <- file.path(dir, "filtered_feature_bc_matrix", "matrix.mtx.gz")
  barcodes_path <- file.path(dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  features_path <- file.path(dir, "filtered_feature_bc_matrix", "features.tsv.gz")
  h5_path <- h5_path %||% file.path(dir, "filtered_feature_bc_matrix", "filtered_feature_bc_matrix.h5")
  positions_path <- file.path(dir, "spatial", "tissue_positions.csv")

  if (!is.null(h5_path) && file.exists(h5_path)) {
    h5_data <- .blade_read_h5_matrix(h5_path)
    matrix <- h5_data$matrix
    barcodes_df <- h5_data$barcodes
    features_df <- h5_data$features
  } else {
    matrix <- Matrix::readMM(gzfile(matrix_path))
    matrix <- Matrix::as(matrix, "dgCMatrix")

    barcodes_df <- polars::read_csv(barcodes_path, separator = "\t", has_header = FALSE)
    features_df <- polars::read_csv(features_path, separator = "\t", has_header = FALSE)
  }
  positions_df <- polars::read_csv(positions_path)

  positions_df <- .blade_standardize_positions(positions_df)

  neighbor_index <- .blade_neighbor_index(positions_df, neighbor_mode, neighbor_radius, k_neighbors)

  list(
    matrix = matrix,
    barcodes = barcodes_df,
    features = features_df,
    positions = positions_df,
    inclusion = .blade_positions_inclusion(positions_df),
    depth = NULL,
    neighbor_index = neighbor_index,
    neighbor_pairs = NULL
  )
}

#' Read Visium HD inputs
#'
#' @param data_dir Optional data directory for auto-discovery.
#' @param h5_path Optional path to filtered_feature_cell_matrix.h5 or filtered_feature_bc_matrix.h5.
#' @param positions_path Optional path to tissue_positions.parquet or tissue_positions.csv.
#' @param neighbor_mode "radius" or "k".
#' @param k_neighbors Number of neighbors for k mode.
#' @param neighbor_radius Optional radius for neighbor search.
#' @return A list with matrix, positions, inclusion mask, and neighbor index.
#' @export
blade_read_visium_hd <- function(
    data_dir = NULL,
    h5_path = NULL,
    positions_path = NULL,
    neighbor_mode = "radius",
    k_neighbors = 8,
    neighbor_radius = NULL) {
  h5_path <- h5_path %||% .blade_auto_find_h5(data_dir)
  if (is.null(h5_path)) {
    stop("HDF5 path not provided and auto-discovery failed.")
  }

  positions_path <- positions_path %||% .blade_auto_find_positions(data_dir, h5_path)
  if (is.null(positions_path)) {
    stop("Positions path not provided and auto-discovery failed.")
  }

  h5_data <- .blade_read_h5_matrix(h5_path)
  matrix <- h5_data$matrix
  barcodes_df <- h5_data$barcodes
  features_df <- h5_data$features

  if (grepl("\\.parquet$", positions_path, ignore.case = TRUE)) {
    positions_df <- polars::read_parquet(positions_path)
  } else {
    positions_df <- polars::read_csv(positions_path)
  }

  positions_df <- .blade_standardize_positions(positions_df)

  neighbor_index <- .blade_neighbor_index(positions_df, neighbor_mode, neighbor_radius, k_neighbors)

  list(
    matrix = matrix,
    barcodes = barcodes_df,
    features = features_df,
    positions = positions_df,
    inclusion = .blade_positions_inclusion(positions_df),
    depth = NULL,
    neighbor_index = neighbor_index,
    neighbor_pairs = NULL
  )
}

#' Remove border spots
#'
#' @param tissue A BLADE tissue list from blade_read_visium() or blade_read_visium_hd().
#' @param distance Border width (for extreme strategy).
#' @param strategy "extreme" or "degree".
#' @param expected_neighbors Expected neighbors for degree strategy.
#' @return Updated tissue list.
#' @export
blade_remove_border <- function(tissue, distance = 0.0, strategy = "extreme", expected_neighbors = 8) {
  border_df <- .blade_get_border(tissue, strategy = strategy, expected_neighbors = expected_neighbors, border_width = distance)
  tissue$inclusion <- .blade_apply_inclusion_zero(tissue$inclusion, border_df$barcode)
  tissue
}

#' Remove edge spots
#'
#' @param tissue A BLADE tissue list.
#' @param distance Edge depth to remove.
#' @param strategy "extreme" or "degree" for border seeding.
#' @param expected_neighbors Expected neighbors for degree strategy.
#' @param seed_mode "zeros", "border", or "both".
#' @param border_width Border width for extreme strategy.
#' @return Updated tissue list.
#' @export
blade_remove_edge <- function(
    tissue,
    distance = 1,
    strategy = "extreme",
    expected_neighbors = 8,
    seed_mode = "zeros",
    border_width = NULL) {
  if (distance <= 0) {
    return(tissue)
  }
  tissue <- .blade_remove_edge_layers(
    tissue,
    distance = distance,
    strategy = strategy,
    expected_neighbors = expected_neighbors,
    seed_mode = seed_mode,
    border_width = border_width
  )
  tissue
}

#' Remove malfunction (outlier) spots
#'
#' @param tissue A BLADE tissue list.
#' @param out_rate Expected outlier rate.
#' @param use_clusters If TRUE, only remove large outlier clusters.
#' @param min_cluster_size Minimum cluster size for clustered removal.
#' @param radius Neighbor radius for clustering.
#' @param radius_scale Scale factor for radius.
#' @return Updated tissue list.
#' @export
blade_remove_malfunction <- function(
    tissue,
    out_rate = 0.05,
    use_clusters = TRUE,
    min_cluster_size = 15,
    radius = NULL,
    radius_scale = 1.0) {
  outliers <- .blade_outlier(.blade_sum_df(tissue), out_rate)
  if (nrow(outliers) == 0) {
    return(tissue)
  }

  if (!use_clusters) {
    tissue$inclusion <- .blade_apply_inclusion_zero(tissue$inclusion, outliers$barcode)
    return(tissue)
  }

  coords <- .blade_positions_df(tissue$positions)
  coords <- coords[coords$barcode %in% outliers$barcode, , drop = FALSE]
  if (nrow(coords) == 0) {
    return(tissue)
  }

  rad <- radius
  if (is.null(rad)) {
    rad <- tissue$neighbor_index$radius
  }
  if (is.null(rad) || !is.finite(rad)) {
    rad <- 1.0
  }
  rad <- rad * radius_scale

  clusters <- .blade_cluster_points(coords$x_mtx, coords$y_mtx, rad)
  to_remove <- unlist(lapply(clusters, function(idx) {
    if (length(idx) >= min_cluster_size) coords$barcode[idx] else character(0)
  }))

  tissue$inclusion <- .blade_apply_inclusion_zero(tissue$inclusion, to_remove)
  tissue
}

#' Iterative edge cleanup loop
#'
#' @param data_dir Visium data directory.
#' @param output_dir Directory to save updated tissue_positions.csv snapshots.
#' @param p_threshold Stop when edge vs interior p-value >= threshold.
#' @param initial_edge_distance Starting edge distance.
#' @param max_iterations Optional maximum iterations.
#' @param remove_border_first Whether to remove border before edge testing.
#' @return A data.frame of edge distances and p-values.
#' @export
blade_run_iterative_edge_cleanup <- function(
    data_dir,
    output_dir,
    p_threshold = 0.05,
    initial_edge_distance = 1,
    max_iterations = NULL,
    remove_border_first = TRUE) {
  tissue <- blade_read_visium(data_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (remove_border_first) {
    tissue <- blade_remove_border(tissue)
    blade_save_positions(tissue, file.path(output_dir, "tissue_positions.border_removed.csv"))
  }
  history <- .blade_iterative_edge_cleanup(
    tissue,
    output_dir = output_dir,
    p_threshold = p_threshold,
    initial_edge_distance = initial_edge_distance,
    max_iterations = max_iterations
  )
  history
}

#' Iterative edge cleanup for Visium HD
#'
#' @param data_dir Optional data directory for auto-discovery.
#' @param h5_path Optional path to a 10x HDF5 matrix file.
#' @param positions_path Optional path to tissue_positions.parquet or CSV.
#' @param output_dir Directory to save tissue_positions snapshots.
#' @param p_threshold Stop when edge vs interior p-value >= threshold.
#' @param initial_edge_distance Starting edge distance.
#' @param max_iterations Optional maximum iterations.
#' @param remove_border_first Whether to remove border before edge testing.
#' @return A data.frame of edge distances and p-values.
#' @export
blade_run_iterative_edge_cleanup_hd <- function(
    data_dir = NULL,
    h5_path = NULL,
    positions_path = NULL,
    output_dir,
    p_threshold = 0.05,
    initial_edge_distance = 1,
    max_iterations = NULL,
    remove_border_first = TRUE) {
  tissue <- blade_read_visium_hd(
    data_dir = data_dir,
    h5_path = h5_path,
    positions_path = positions_path
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (remove_border_first) {
    tissue <- blade_remove_border(tissue)
    blade_save_positions(tissue, file.path(output_dir, "tissue_positions.border_removed.csv"))
  }

  .blade_iterative_edge_cleanup(
    tissue,
    output_dir = output_dir,
    p_threshold = p_threshold,
    initial_edge_distance = initial_edge_distance,
    max_iterations = max_iterations
  )
}

#' Save updated tissue positions
#'
#' @param tissue A BLADE tissue list.
#' @param out_path Output CSV path.
#' @return Invisibly returns the path.
#' @export
blade_save_positions <- function(tissue, out_path) {
  pos_pl <- tissue$positions
  inc_pl <- tissue$inclusion
  merged <- pos_pl$join(inc_pl, on = "barcode", how = "left")
  out_df <- merged$to_data_frame()
  if ("in_tissue_right" %in% names(out_df)) {
    out_df$in_tissue <- ifelse(is.na(out_df$in_tissue_right), out_df$in_tissue, out_df$in_tissue_right)
    out_df$in_tissue_right <- NULL
  }
  utils::write.csv(out_df, out_path, row.names = FALSE)
  invisible(out_path)
}

# Helpers

`%||%` <- function(x, y) if (is.null(x)) y else x

.blade_barcodes_vec <- function(barcodes_df) {
  if ("barcode" %in% barcodes_df$columns) {
    return(as.character(barcodes_df$get_column("barcode")$to_list()))
  }
  as.character(barcodes_df$get_column(barcodes_df$columns[[1]])$to_list())
}

.blade_positions_inclusion <- function(positions_df) {
  polars::pl$DataFrame(list(
    barcode = positions_df$get_column("barcode")$to_list(),
    in_tissue = as.integer(positions_df$get_column("in_tissue")$to_list())
  ))
}

.blade_standardize_positions <- function(positions_df) {
  cols <- positions_df$columns
  id_col <- NULL
  for (cand in c("barcode", "cell_id", "bin_id", "barcode_id")) {
    if (cand %in% cols) {
      id_col <- cand
      break
    }
  }
  if (is.null(id_col)) {
    stop("Positions file must contain a barcode/cell identifier column.")
  }

  coord_pairs <- list(
    c("x_mtx", "y_mtx"),
    c("x", "y"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres")
  )
  x_col <- y_col <- NULL
  for (pair in coord_pairs) {
    if (all(pair %in% cols)) {
      x_col <- pair[1]
      y_col <- pair[2]
      break
    }
  }
  if (is.null(x_col)) {
    stop("Positions file must contain coordinate columns (e.g. x_mtx/y_mtx).")
  }

  if (!("in_tissue" %in% cols)) {
    positions_df <- positions_df$with_columns(polars::pl$lit(1L)$alias("in_tissue"))
  }

  positions_df$select(c(
    polars::pl$col(id_col)$cast(polars::pl$Utf8)$alias("barcode"),
    polars::pl$col("in_tissue")$cast(polars::pl$Int64)$alias("in_tissue"),
    polars::pl$col(x_col)$cast(polars::pl$Float64)$alias("x_mtx"),
    polars::pl$col(y_col)$cast(polars::pl$Float64)$alias("y_mtx"),
    polars::pl$col(x_col)$cast(polars::pl$Float64)$alias("x_coor"),
    polars::pl$col(y_col)$cast(polars::pl$Float64)$alias("y_coor")
  ))
}

.blade_auto_find_h5 <- function(data_dir) {
  if (is.null(data_dir)) return(NULL)
  candidates <- c(
    file.path(data_dir, "filtered_feature_cell_matrix.h5"),
    file.path(data_dir, "filtered_feature_bc_matrix.h5")
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  NULL
}

.blade_auto_find_positions <- function(data_dir, h5_path) {
  root <- if (!is.null(data_dir)) data_dir else dirname(h5_path)
  candidates <- c(
    file.path(root, "spatial", "tissue_positions.parquet"),
    file.path(root, "spatial", "tissue_positions.csv")
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  NULL
}

.blade_read_h5_matrix <- function(h5_path) {
  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(h5$close_all())
  grp <- h5[["matrix"]]
  data <- grp[["data"]][]
  indices <- grp[["indices"]][]
  indptr <- grp[["indptr"]][]
  shape <- as.integer(grp[["shape"]][])
  barcodes <- grp[["barcodes"]][]

  feat_grp <- grp[["features"]]
  feature_id <- feat_grp[["id"]][]
  feature_name <- feat_grp[["name"]][]
  feature_type <- if ("feature_type" %in% names(feat_grp)) feat_grp[["feature_type"]][] else NULL

  matrix <- Matrix::sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = data,
    dims = shape,
    repr = "C"
  )

  barcodes_df <- polars::pl$DataFrame(list(barcode = as.character(barcodes)))
  features_list <- list(feature_id = as.character(feature_id), feature_name = as.character(feature_name))
  if (!is.null(feature_type)) {
    features_list$feature_type <- as.character(feature_type)
  }
  features_df <- polars::pl$DataFrame(features_list)
  list(matrix = matrix, barcodes = barcodes_df, features = features_df)
}

.blade_neighbor_index <- function(positions_df, mode, radius, k) {
  radius <- radius
  if (is.null(radius)) {
    radius <- .blade_estimate_grid_radius(positions_df)
  }
  list(positions = positions_df, mode = mode, radius = radius, k = k, cache = new.env(parent = emptyenv()))
}

.blade_estimate_grid_radius <- function(pos_pl) {
  median_step <- function(values) {
    uniq <- sort(unique(values))
    if (length(uniq) < 2) return(NA_real_)
    diffs <- diff(uniq)
    diffs <- diffs[diffs > 0]
    if (length(diffs) == 0) return(NA_real_)
    median(diffs)
  }
  x_vals <- pos_pl$get_column("x_mtx")$to_list()
  y_vals <- pos_pl$get_column("y_mtx")$to_list()
  dx <- median_step(x_vals)
  dy <- median_step(y_vals)
  steps <- c(dx, dy)
  steps <- steps[is.finite(steps)]
  if (length(steps) == 0) return(NA_real_)
  max(steps) * 2.1
}

.blade_neighbors_by_barcode <- function(index, barcode, radius = NULL, k = NULL, include_center = FALSE) {
  rad <- radius
  if (is.null(rad)) {
    if (!is.null(k)) {
      rad <- index$radius
    } else if (index$mode == "k") {
      k <- index$k
      rad <- index$radius
    } else {
      rad <- index$radius
    }
  }
  if (is.null(rad) || !is.finite(rad)) {
    return(character())
  }

  key <- paste(rad, k, include_center, sep = "|")
  if (exists(key, envir = index$cache, inherits = FALSE)) {
    cache_val <- get(key, envir = index$cache, inherits = FALSE)
    return(cache_val[[barcode]] %||% character())
  }

  pairs <- .blade_neighbor_pairs(index$positions, rad)
  if (!include_center) {
    pairs <- pairs$filter(polars::pl$col("neighbor") != polars::pl$col("barcode"))
  }
  pairs <- pairs$sort(c("barcode", "distance"))

  if (!is.null(k)) {
    grouped <- pairs$group_by("barcode")$agg(polars::pl$col("neighbor")$head(k))
  } else {
    grouped <- pairs$group_by("barcode")$agg(polars::pl$col("neighbor"))
  }
  grouped_df <- grouped$to_data_frame()
  split_pairs <- setNames(grouped_df$neighbor, grouped_df$barcode)

  assign(key, split_pairs, envir = index$cache)
  split_pairs[[barcode]] %||% character()
}

.blade_neighbor_pairs <- function(pos_pl, radius) {
  if (radius <= 0) {
    return(polars::pl$DataFrame(list(barcode = character(), neighbor = character(), distance = numeric())))
  }
  cell_size <- radius
  left <- pos_pl$with_columns(
    (polars::pl$col("x_mtx") / cell_size)$floor()$cast(polars::pl$Int64)$alias("cell_x"),
    (polars::pl$col("y_mtx") / cell_size)$floor()$cast(polars::pl$Int64)$alias("cell_y")
  )
  offsets <- polars::pl$DataFrame(list(
    dx = c(-1, 0, 1, -1, 0, 1, -1, 0, 1),
    dy = c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
  ))
  candidates <- left$join(offsets, how = "cross")$with_columns(
    (polars::pl$col("cell_x") + polars::pl$col("dx"))$alias("cell_x_n"),
    (polars::pl$col("cell_y") + polars::pl$col("dy"))$alias("cell_y_n")
  )
  right <- left$select(
    polars::pl$col("barcode")$alias("neighbor"),
    polars::pl$col("x_mtx")$alias("nx"),
    polars::pl$col("y_mtx")$alias("ny"),
    polars::pl$col("cell_x"),
    polars::pl$col("cell_y")
  )
  pairs <- candidates$join(right, left_on = c("cell_x_n", "cell_y_n"), right_on = c("cell_x", "cell_y"), how = "inner")$
    with_columns(
      ((polars::pl$col("x_mtx") - polars::pl$col("nx"))^2 + (polars::pl$col("y_mtx") - polars::pl$col("ny"))^2)$sqrt()$alias("distance")
    )$
    filter(polars::pl$col("distance") <= radius)$
    select(c("barcode", "neighbor", "distance"))
  pairs
}

.blade_sum_df <- function(tissue) {
  counts <- Matrix::colSums(tissue$matrix)
  barcodes <- .blade_barcodes_vec(tissue$barcodes)
  pos_pl <- tissue$positions$select(c("barcode", "x_mtx", "y_mtx"))
  counts_pl <- polars::pl$DataFrame(list(barcode = barcodes, gene_count = as.numeric(counts)))
  merged <- pos_pl$join(counts_pl, on = "barcode", how = "left")
  merged$to_data_frame()
}

.blade_outlier <- function(df, out_rate = 0.05) {
  out_rate_epr <- 0
  probs <- c(0.02, 0.98)
  data_outlier <- df
  counts <- df$gene_count
  skewness <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 3) return(0)
    m <- mean(x)
    s <- stats::sd(x)
    if (s == 0) return(0)
    mean(((x - m) / s)^3)
  }
  while (abs(out_rate_epr - out_rate) > 0.01) {
    if (out_rate_epr - out_rate > 0) {
      probs <- c(probs[1] - 0.0005, probs[2] + 0.0005)
    }
    if (out_rate_epr - out_rate < 0) {
      probs <- c(probs[1] + 0.0005, probs[2] - 0.0005)
    }
    quart_1 <- stats::quantile(counts, probs[1])
    quart_2 <- stats::quantile(counts, probs[2])
    iqr <- stats::quantile(counts, 0.75) - stats::quantile(counts, 0.25)
    medcp <- skewness(counts)
    lower <- quart_1 - 1.5 * exp(-4 * medcp) * iqr
    upper <- quart_2 + 1.5 * exp(3 * medcp) * iqr
    mask <- !(counts > lower & counts < upper)
    data_outlier <- df[mask, , drop = FALSE]
    out_rate_epr <- nrow(data_outlier) / nrow(df)
  }
  data_outlier
}

.blade_get_border <- function(tissue, expected_neighbors = 8, strategy = "extreme", border_width = 0.0) {
  pos_pl <- tissue$positions
  inclusion_pl <- tissue$inclusion
  included <- inclusion_pl$filter(polars::pl$col("in_tissue") == 1)$select("barcode")
  covered <- pos_pl$join(included, on = "barcode", how = "inner")
  if (covered$height == 0) return(data.frame(barcode = character()))

  if (strategy == "extreme") {
    xmin <- covered$select(polars::pl$col("x_mtx")$min())$to_list()[[1]][[1]]
    xmax <- covered$select(polars::pl$col("x_mtx")$max())$to_list()[[1]][[1]]
    ymin <- covered$select(polars::pl$col("y_mtx")$min())$to_list()[[1]][[1]]
    ymax <- covered$select(polars::pl$col("y_mtx")$max())$to_list()[[1]][[1]]
    bw <- as.numeric(border_width)
    keep <- covered$filter(
      (polars::pl$col("x_mtx") - xmin)$abs() <= bw |
        (polars::pl$col("x_mtx") - xmax)$abs() <= bw |
        (polars::pl$col("y_mtx") - ymin)$abs() <= bw |
        (polars::pl$col("y_mtx") - ymax)$abs() <= bw
    )
    return(keep$to_data_frame())
  }

  covered_df <- covered$select("barcode")$to_data_frame()
  bars <- character()
  for (i in seq_len(nrow(covered_df))) {
    b <- covered_df$barcode[i]
    ngh <- .blade_neighbors_by_barcode(tissue$neighbor_index, b)
    if (length(ngh) < expected_neighbors) {
      bars <- c(bars, b)
    }
  }
  covered$filter(polars::pl$col("barcode")$is_in(bars))$to_data_frame()
}

.blade_spot_depth <- function(tissue) {
  inclusion_pl <- tissue$inclusion
  zero_barcodes <- inclusion_pl$filter(polars::pl$col("in_tissue") == 0)$get_column("barcode")$to_list()
  covered_barcodes <- inclusion_pl$filter(polars::pl$col("in_tissue") == 1)$get_column("barcode")$to_list()

  frontier <- character()
  for (b in zero_barcodes) {
    frontier <- unique(c(frontier, .blade_neighbors_by_barcode(tissue$neighbor_index, b)))
  }
  border_barcodes <- .blade_get_border(tissue)$barcode
  frontier <- unique(c(frontier, border_barcodes))
  frontier <- intersect(frontier, covered_barcodes)

  depth_map <- setNames(rep(0L, length(covered_barcodes)), covered_barcodes)
  depth_map[frontier] <- 1L
  visited <- frontier
  level <- 2L

  while (length(setdiff(covered_barcodes, visited)) > 0) {
    next_frontier <- character()
    for (b in frontier) {
      next_frontier <- unique(c(next_frontier, .blade_neighbors_by_barcode(tissue$neighbor_index, b)))
    }
    next_frontier <- intersect(next_frontier, covered_barcodes)
    next_frontier <- setdiff(next_frontier, visited)
    if (length(next_frontier) == 0) break
    depth_map[next_frontier] <- level
    visited <- c(visited, next_frontier)
    frontier <- next_frontier
    level <- level + 1L
  }

  data.frame(barcode = names(depth_map), spot_depth = unname(depth_map))
}

.blade_apply_inclusion_zero <- function(inclusion_df, barcodes) {
  if (length(barcodes) == 0) return(inclusion_df)
  inclusion_df$with_columns(
    polars::pl$when(polars::pl$col("barcode")$is_in(barcodes))$
      then(0L)$
      otherwise(polars::pl$col("in_tissue"))$
      alias("in_tissue")
  )
}

.blade_remove_edge_layers <- function(tissue, distance, strategy, expected_neighbors, seed_mode, border_width) {
  inclusion_pl <- tissue$inclusion
  covered_barcodes <- inclusion_pl$filter(polars::pl$col("in_tissue") == 1)$get_column("barcode")$to_list()

  frontier <- character()
  if (seed_mode %in% c("zeros", "both")) {
    zero_barcodes <- inclusion_pl$filter(polars::pl$col("in_tissue") == 0)$get_column("barcode")$to_list()
    margin <- tissue$neighbor_index$radius
    for (b in zero_barcodes) {
      frontier <- unique(c(frontier, .blade_neighbors_by_barcode(tissue$neighbor_index, b)))
    }
    if (!is.null(margin) && is.finite(margin) && length(frontier) > 0 && length(covered_barcodes) > 0) {
      covered_pl <- tissue$positions$filter(polars::pl$col("barcode")$is_in(covered_barcodes))
      xmin <- covered_pl$select(polars::pl$col("x_mtx")$min())$to_list()[[1]][[1]]
      xmax <- covered_pl$select(polars::pl$col("x_mtx")$max())$to_list()[[1]][[1]]
      ymin <- covered_pl$select(polars::pl$col("y_mtx")$min())$to_list()[[1]][[1]]
      ymax <- covered_pl$select(polars::pl$col("y_mtx")$max())$to_list()[[1]][[1]]
      band <- margin * 2.0
      filtered <- tissue$positions$
        filter(polars::pl$col("barcode")$is_in(frontier))$
        filter(
          polars::pl$col("x_mtx") >= (xmin - band) &
            polars::pl$col("x_mtx") <= (xmax + band) &
            polars::pl$col("y_mtx") >= (ymin - band) &
            polars::pl$col("y_mtx") <= (ymax + band)
        )$
        get_column("barcode")$to_list()
      frontier <- filtered
    }
  }

  if (seed_mode %in% c("border", "both")) {
    bw <- border_width
    if (is.null(bw)) {
      bw <- tissue$neighbor_index$radius
    }
    border_df <- .blade_get_border(tissue, strategy = strategy, expected_neighbors = expected_neighbors, border_width = bw)
    frontier <- unique(c(frontier, border_df$barcode))
  }

  frontier <- intersect(frontier, covered_barcodes)
  visited <- frontier
  current_frontier <- frontier
  if (distance > 1) {
    for (i in seq_len(distance - 1)) {
      next_frontier <- character()
      for (b in current_frontier) {
        next_frontier <- unique(c(next_frontier, .blade_neighbors_by_barcode(tissue$neighbor_index, b)))
      }
      next_frontier <- intersect(next_frontier, covered_barcodes)
      next_frontier <- setdiff(next_frontier, visited)
      if (length(next_frontier) == 0) break
      visited <- unique(c(visited, next_frontier))
      current_frontier <- next_frontier
    }
  }

  tissue$inclusion <- .blade_apply_inclusion_zero(tissue$inclusion, visited)
  tissue
}

.blade_iterative_edge_cleanup <- function(tissue, output_dir, p_threshold, initial_edge_distance, max_iterations) {
  barcodes <- .blade_barcodes_vec(tissue$barcodes)
  gene_counts <- Matrix::colSums(tissue$matrix)
  depth <- .blade_spot_depth(tissue)
  depth_map <- setNames(depth$spot_depth, depth$barcode)
  depth_values <- as.integer(depth_map[barcodes])
  inclusion_map <- setNames(tissue$inclusion$get_column("in_tissue")$to_list(), tissue$inclusion$get_column("barcode")$to_list())

  history <- data.frame(edge_distance = integer(), p_value = numeric())
  current_distance <- max(1L, as.integer(initial_edge_distance))
  iteration <- 0

  while (current_distance <= max(depth_values, na.rm = TRUE)) {
    mask <- inclusion_map[barcodes] == 1
    edge_mask <- depth_values == current_distance & mask
    interior_mask <- depth_values > current_distance & mask
    if (!any(edge_mask) || !any(interior_mask)) {
      break
    }

    p_val <- stats::t.test(gene_counts[edge_mask], gene_counts[interior_mask])$p.value
    history <- rbind(history, data.frame(edge_distance = current_distance, p_value = p_val))

    if (p_val < p_threshold) {
      tissue <- .blade_remove_edge_layers(tissue, distance = current_distance, strategy = "extreme", expected_neighbors = 8, seed_mode = "both", border_width = NULL)
      inclusion_map <- setNames(tissue$inclusion$get_column("in_tissue")$to_list(), tissue$inclusion$get_column("barcode")$to_list())
      iteration <- iteration + 1
      blade_save_positions(tissue, file.path(output_dir, sprintf("tissue_positions.iter%d.csv", iteration)))
      current_distance <- current_distance + 1
      if (!is.null(max_iterations) && iteration >= max_iterations) {
        break
      }
    } else {
      break
    }
  }

  history
}

.blade_cluster_points <- function(x, y, radius) {
  n <- length(x)
  visited <- rep(FALSE, n)
  clusters <- list()

  for (i in seq_len(n)) {
    if (visited[i]) next
    queue <- i
    visited[i] <- TRUE
    cluster <- i
    while (length(queue) > 0) {
      idx <- queue[1]
      queue <- queue[-1]
      dists <- sqrt((x - x[idx])^2 + (y - y[idx])^2)
      neighbors <- which(dists <= radius)
      new_neighbors <- neighbors[!visited[neighbors]]
      if (length(new_neighbors) > 0) {
        visited[new_neighbors] <- TRUE
        queue <- c(queue, new_neighbors)
        cluster <- c(cluster, new_neighbors)
      }
    }
    clusters[[length(clusters) + 1]] <- cluster
  }

  clusters
}
