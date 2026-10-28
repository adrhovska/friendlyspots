#' Spatial similarity between spatial transcriptomics samples
#'
#' @param samples Named list of data frames, one per sample, with columns
#'   `x`, `y` (positions in micrometres) and `label`.
#' @param normalize Divide each sample's connection counts by their total.
#' @param pools Optional named list mapping pool names to sample names;
#'   connection counts are summed within each pool.
#' @param radius Neighbourhood radius in micrometres. The default
#'   includes the six immediate Visium neighbours (~100 um apart, allowing for
#'   position jitter) but not the next ring at ~173 um.
#' @return Symmetric similarity matrix with values in \[0, 1\].
#' @examples
#' samples <- list(
#'   a = data.frame(x = c(0, 100, 50), y = c(0, 0, 86.6),
#'                  label = c("x", "y", "x")),
#'   b = data.frame(x = c(0, 100, 200), y = c(0, 0, 0),
#'                  label = c("x", "x", "y")),
#'   c = data.frame(x = c(0, 100), y = c(0, 0),
#'                  label = c("y", "y"))
#' )
#' spatial_similarity(samples)
#' spatial_similarity(samples, pools = list(ab = c("a", "b"), c = "c"))
#' @importFrom stats dist
#' @export
spatial_similarity <- function(samples, normalize = TRUE, pools = NULL,
                               radius = 105) {
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
  counts <- do.call(rbind, lapply(samples, connection_counts, labels = labels,
                                   radius = radius))
  if (!is.null(pools)) {
    counts <- do.call(rbind, lapply(pools, function(p) {
      colSums(counts[p, , drop = FALSE])
    }))
  }
  if (normalize) counts <- counts / rowSums(counts)
  d <- as.matrix(dist(counts))
  1 - d / max(d)
}

connection_counts <- function(sample, labels, radius) {
  label <- factor(sample$label, levels = labels)
  nb <- as.matrix(dist(sample[, c("x", "y")])) <= radius
  diag(nb) <- FALSE
  onehot <- outer(as.integer(label), seq_along(labels), "==") * 1
  counts <- t(onehot) %*% nb %*% onehot
  as.vector(t(counts))
}