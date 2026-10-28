import warnings

import numpy as np
import pandas as pd

# Hexagonal Visium neighbors in array coordinates.
HEX_OFFSETS = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]


def spatial_similarity(samples, normalize=True, pools=None):
    """Spatial similarity between spatial transcriptomics samples.

    samples: dict of DataFrames with columns array_row, array_col (hexagonal
        grid array coordinates) and label.
    normalize: divide each sample's connection counts by their total.
    pools: optional dict mapping pool names to lists of sample names;
        connection counts are summed within each pool.
    Returns a symmetric DataFrame of similarities in [0, 1].
    """
    if any(s["label"].isna().any() for s in samples.values()):
        # e.g. pd.read_csv() reads a label literally called "NA" as missing
        warnings.warn(
            "Some labels are 'NA'. These will be treated as a label called "
            "'NA' and used in the comparison, rather than being treated as "
            "the absence of a label and dropped."
        )
    samples = {name: s.assign(label=s["label"].astype(str).where(s["label"].notna(), "NA")) for name, s in samples.items()}
    # Sorted as set order varies between runs, changing floating-point results
    labels = sorted(set().union(*(s["label"].astype(str) for s in samples.values())))
    counts = pd.DataFrame({name: _connection_counts(s, labels) for name, s in samples.items()}).T
    if pools is not None:
        counts = pd.DataFrame({name: counts.loc[members].sum() for name, members in pools.items()}).T
    if normalize:
        counts = counts.div(counts.sum(axis=1), axis=0)
    x = counts.to_numpy(dtype=float)
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
    return pd.DataFrame(1 - d / d.max(), index=counts.index, columns=counts.index)


def _connection_counts(sample, labels):
    codes = pd.Categorical(sample["label"].astype(str), categories=labels).codes
    coords = list(zip(sample["array_row"], sample["array_col"]))
    index = {rc: i for i, rc in enumerate(coords)}
    counts = np.zeros((len(labels), len(labels)))
    for i, (r, c) in enumerate(coords):
        for dr, dc in HEX_OFFSETS:
            j = index.get((r + dr, c + dc))
            if j is not None:
                counts[codes[i], codes[j]] += 1
    return counts.ravel()
