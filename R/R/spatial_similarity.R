#' Spatial similarity between spatial transcriptomics samples
#'
#' @param samples Named list of data frames, one per sample, with columns
#'   `array_row`, `array_col` (hexagonal grid array coordinates) and `label`.
#' @param normalize Divide each sample's connection counts by their total.
#' @param pools Optional named list mapping pool names to sample names;
#'   connection counts are summed within each pool.
#' @return Symmetric similarity matrix with values in \[0, 1\].
#' @examples
#' samples <- list(
#'   a = data.frame(array_row = c(0, 0, 1), array_col = c(0, 2, 1),
#'                  label = c("x", "y", "x")),
#'   b = data.frame(array_row = c(0, 0, 0), array_col = c(0, 2, 4),
#'                  label = c("x", "x", "y")),
#'   c = data.frame(array_row = c(0, 0), array_col = c(0, 2),
#'                  label = c("y", "y"))
#' )
#' spatial_similarity(samples)
#' spatial_similarity(samples, pools = list(ab = c("a", "b"), c = "c"))
#' @importFrom stats dist
#' @export
spatial_similarity <- function(samples, normalize = TRUE, pools = NULL) {
  if (anyNA(unlist(lapply(samples, function(s) as.character(s$label))))) {
    # e.g. read.csv() reads a label literally called "NA" as missing
    warning("Some labels are 'NA'. These will be treated as a label called ",
            "'NA' and used in the comparison, rather than being treated as ",
            "the absence of a label and dropped.")
  }
  samples <- lapply(samples, function(s) {
    s$label <- ifelse(is.na(s$label), "NA", as.character(s$label))
    s
  })
  # Sorted so floating-point results don't depend on sample order
  labels <- sort(unique(unlist(lapply(samples, function(s) s$label))))
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
