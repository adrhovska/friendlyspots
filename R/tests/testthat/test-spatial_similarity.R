samples <- list(
  A = data.frame(array_row = c(0, 0, 1), array_col = c(0, 2, 1), label = c("a", "b", "a")),
  B = data.frame(array_row = c(0, 0, 0, 2), array_col = c(0, 2, 4, 0), label = c("a", "a", "b", "b")),
  C = data.frame(array_row = c(0, 0), array_col = c(0, 2), label = c("b", "b"))
)

expect_sim <- function(sim, ab, ac, bc) {
  n <- rownames(sim)
  expected <- matrix(c(1, ab, ac, ab, 1, bc, ac, bc, 1), 3, dimnames = list(n, n))
  expect_equal(sim, expected)
}

test_that("raw counts give 1 - d / max(d)", {
  sim <- spatial_similarity(samples, normalize = FALSE)
  expect_sim(sim, 1 - sqrt(2) / 4, 0, 1 - sqrt(10) / 4)
})

test_that("normalize divides counts by their total", {
  dmax <- sqrt(11 / 8)
  sim <- spatial_similarity(samples)
  expect_sim(sim, 1 - (sqrt(6) / 12) / dmax, 1 - sqrt(4 / 3) / dmax, 0)
})

test_that("pools sum counts of their samples", {
  pools <- list(p1 = c("A", "B"), p2 = "C", p3 = "A")
  sim <- spatial_similarity(samples, normalize = FALSE, pools = pools)
  expect_equal(rownames(sim), c("p1", "p2", "p3"))
  expect_sim(sim, 0, 1 - sqrt(6) / sqrt(38), 1 - 4 / sqrt(38))
})
