import numpy as np
import pandas as pd


def spatial_similarity(samples, normalize=True, pools=None, radius=105):
    # 105 because the nearest neighbours in Visium ~ 100
    # and this margin makes sure all are included. The next layer of neighbours is at 173 um.
    """Spatial similarity between spatial transcriptomics samples.

    samples: dict of DataFrames with columns x, y (positions in micrometres) and label.
    normalize: divide each sample's connection counts by their total.
    pools: optional dict mapping pool names to lists of sample names;
        connection counts are summed within each pool.
    radius: neighbourhood radius in micrometres. The default includes the six
        immediate Visium neighbours (100 um apart), allows for flexibility in technology
        used and for altering radius value.
    Returns a symmetric DataFrame of similarities in [0, 1].
    """
    labels = sorted(set().union(*(s["label"].astype(str) for s in samples.values())))
    counts = pd.DataFrame(
        {name: _connection_counts(s, labels, radius) for name, s in samples.items()}
    ).T
    if pools is not None:
        counts = pd.DataFrame(
            {name: counts.loc[members].sum() for name, members in pools.items()}
        ).T
    if normalize:
        counts = counts.div(counts.sum(axis=1), axis=0)
    x = counts.to_numpy(dtype=float)
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
    return pd.DataFrame(1 - d / d.max(), index=counts.index, columns=counts.index)


def _connection_counts(sample, labels, radius):
    codes = pd.Categorical(sample["label"].astype(str), categories=labels).codes
    xy = sample[["x", "y"]].to_numpy(dtype=float)
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    neighbours = d <= radius
    np.fill_diagonal(neighbours, False)
    onehot = np.eye(len(labels))[codes]
    return (onehot.T @ neighbours @ onehot).ravel()
