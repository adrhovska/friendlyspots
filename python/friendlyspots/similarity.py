import warnings

import numpy as np
import pandas as pd

# Hexagonal Visium neighbors in array coordinates.
HEX_OFFSETS = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]


def spatial_similarity(
    samples, sample_normalize=True, label_normalize=False, pools=None, only_transitions=False, normalize=None
):
    """Spatial similarity between spatial transcriptomics samples.

    samples: dict of DataFrames with columns array_row, array_col (hexagonal
        grid array coordinates) and label.
    sample_normalize: divide each sample's connection counts by their total.
    label_normalize: control for the size of each label: add one connection
        of every label to itself, divide each count by the square root of the
        two labels' connections to themselves (so each label's connections to
        itself become 1), then divide each label's counts by their sum (its
        degree), so each label's counts become the proportions of its
        neighbors' labels, as in one step of a random walk. Applied after
        pooling. Use instead of sample_normalize: set sample_normalize=False.
    pools: optional dict mapping pool names to lists of sample names;
        connection counts are summed within each pool.
    only_transitions: if True, only count connections between spots with
        distinct labels (i.e. borders); does not count connections between
        spots with the same label.
    normalize: deprecated; use sample_normalize.
    Returns a symmetric DataFrame of similarities in [0, 1], or NaN for
    samples with no neighboring spots when sample_normalize=True.
    """
    if normalize is not None:
        warnings.warn("normalize is deprecated; use sample_normalize.", DeprecationWarning, stacklevel=2)
        sample_normalize = normalize
    if sample_normalize and label_normalize:
        raise ValueError(
            "Use either sample_normalize or label_normalize, not both: set sample_normalize=False to use label_normalize."
        )
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
    counts = pd.DataFrame({name: _connection_counts(s, labels, only_transitions) for name, s in samples.items()}).T
    if pools is not None:
        counts = pd.DataFrame({name: counts.loc[members].sum() for name, members in pools.items()}).T
    if label_normalize:
        x = counts.to_numpy(dtype=float).reshape(-1, len(labels), len(labels)) + np.eye(len(labels))
        d = np.sqrt(np.diagonal(x, axis1=1, axis2=2))
        x = x / (d[:, :, None] * d[:, None, :])
        x = x / x.sum(axis=2, keepdims=True)
        counts = pd.DataFrame(x.reshape(len(counts), -1), index=counts.index, columns=counts.columns)
    if sample_normalize:
        totals = counts.sum(axis=1)
        if (totals == 0).any():
            warnings.warn(
                "Similarities are NA for samples with no neighboring spots"
                + (" with different labels (no transitions)" if only_transitions else "")
                + ": "
                + ", ".join(map(str, counts.index[totals == 0]))
            )
        counts = counts.div(totals, axis=0)
    x = counts.to_numpy(dtype=float)
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
    # Identical samples give a max distance of 0, so don't divide by it
    d_max = np.nanmax(d, initial=0)
    if d_max > 0:
        d = d / d_max
    return pd.DataFrame(1 - d, index=counts.index, columns=counts.index)


def _connection_counts(sample, labels, only_transitions=False):
    codes = pd.Categorical(sample["label"].astype(str), categories=labels).codes
    coords = list(zip(sample["array_row"], sample["array_col"]))
    index = {rc: i for i, rc in enumerate(coords)}
    counts = np.zeros((len(labels), len(labels)))
    for i, (r, c) in enumerate(coords):
        for dr, dc in HEX_OFFSETS:
            j = index.get((r + dr, c + dc))
            if j is not None:
                counts[codes[i], codes[j]] += 1
    if only_transitions:
        np.fill_diagonal(counts, 0)
    return counts.ravel()
