# friendlyspots

Compute similarity between spatial transcriptomics samples based on which
labels (e.g. cell types) border one another.

## Installation

```r
remotes::install_github("george-hall-ucl/friendlyspots", subdir = "R")
```

Requires only base R. Also available as a
[Python package](https://github.com/george-hall-ucl/friendlyspots/tree/main/python).

## Usage

`friendlyspots` is only compatible with original (non-HD) Visium.

```r
library(friendlyspots)

# Create one data frame per sample, with columns for array_row, array_col,
# label.
#   e.g. starting from a sample's Seurat object `cells` (requires a VisiumV1
#   image, the default before Seurat 5.1):
#   coords <- cells@images[[1]]@coordinates
#   df1 <- data.frame(array_row = coords$row, array_col = coords$col,
#                     label = cells$cell_type[rownames(coords)])
samples <- list(sample_1 = df1, sample_2 = df2, sample_3 = df3)

sim <- spatial_similarity(samples)
```

Returns a sample × sample matrix with normalized similarity scores: 1 =
identical, 0 = least similar pair.

### Options

- `normalize = FALSE`: compare raw connection counts instead of proportions.
  This will make samples with different numbers of spots be considered
  spatially different even if their underlying spatial features are similar.
- `pools = list(group_1 = c("sample_1", "sample_2"), group_2 = "sample_3")`:
  sum connection counts within groups before comparing.

## Citation

The algorithm and original implementation of `friendlyspots` are based on the
methodology in:

- Durkin _et al_. "Functional integration of an autologous engineered esophagus
  in a large-animal model." _Nature Biotechnology_ (2026).
- Budinger _et al_. "An in vivo and in vitro spatiotemporal profile of human
  midbrain development." _Nature Communications_ (2026).

## License

Copyright (C) 2026 University College London. Licensed under GPL-3.0-or-later.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
