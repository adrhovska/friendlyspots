#' Spatial similarity between spatial transcriptomics samples
#'
#' @param samples Named list of data frames, one per sample, with columns
#'   `array_row`, `array_col` (hexagonal grid array coordinates) and `label`.
#' @param normalize Divide each sample's connection counts by their total.
#' @param pools Optional named list mapping pool names to sample names;
#'   connection counts are summed within each pool.
#' @return Symmetric similarity matrix with values in [0, 1].
#' @export
spatial_similarity <- function(samples, normalize = TRUE, pools = NULL) {
  labels <- sort(unique(unlist(lapply(samples, function(s) as.character(s$label)))))
  counts <- do.call(rbind, lapply(samples, connection_counts, labels = labels))
  if (!is.null(pools)) {
    counts <- do.call(rbind, lapply(pools, function(p) colSums(counts[p, , drop = FALSE])))
  }
  if (normalize) counts <- counts / rowSums(counts)
  d <- as.matrix(dist(counts))
  1 - d / max(d)
}

# Hexagonal Visium neighbors in array coordinates.
HEX_OFFSETS <- list(c(0, -2), c(0, 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))

connection_counts <- function(sample, labels) {
  label <- factor(sample$label, levels = labels)
  key <- paste(sample$array_row, sample$array_col)
  counts <- matrix(0, length(labels), length(labels))
  for (o in HEX_OFFSETS) {
    nbr <- match(paste(sample$array_row + o[1], sample$array_col + o[2]), key)
    hit <- !is.na(nbr)
    counts <- counts + unclass(table(label[hit], label[nbr[hit]]))
  }
  as.vector(t(counts))
}
