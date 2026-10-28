import io

import numpy as np
import pandas as pd
import pytest

from friendlyspots import spatial_similarity

SAMPLES = {
    "A": pd.DataFrame({"array_row": [0, 0, 1], "array_col": [0, 2, 1], "label": ["a", "b", "a"]}),
    "B": pd.DataFrame({"array_row": [0, 0, 0, 2], "array_col": [0, 2, 4, 0], "label": ["a", "a", "b", "b"]}),
    "C": pd.DataFrame({"array_row": [0, 0], "array_col": [0, 2], "label": ["b", "b"]}),
}


def assert_sim(sim, ab, ac, bc):
    expected = np.array([[1, ab, ac], [ab, 1, bc], [ac, bc, 1]])
    np.testing.assert_allclose(sim.to_numpy(), expected)
    assert list(sim.columns) == list(sim.index)


def test_raw_counts():
    sim = spatial_similarity(SAMPLES, normalize=False)
    assert_sim(sim, 1 - np.sqrt(2) / 4, 0, 1 - np.sqrt(10) / 4)
    # pd.read_csv() turns labels written as NA, nan, None, etc. into NaN
    csv = pd.read_csv(io.StringIO(
        "sample,array_row,array_col,label\n"
        "A,0,0,a\nA,0,2,NA\nA,1,1,a\n"
        "B,0,0,a\nB,0,2,a\nB,0,4,nan\nB,2,0,nan\n"
        "C,0,0,None\nC,0,2,None\n"
    ))
    from_file = {name: s.drop(columns="sample") for name, s in csv.groupby("sample")}
    with pytest.warns(UserWarning, match="treated as a label called 'NA'"):
        pd.testing.assert_frame_equal(spatial_similarity(from_file, normalize=False), sim)


def test_normalize():
    dmax = np.sqrt(11 / 8)
    sim = spatial_similarity(SAMPLES)
    assert_sim(sim, 1 - (np.sqrt(6) / 12) / dmax, 1 - np.sqrt(4 / 3) / dmax, 0)


def test_pools():
    pools = {"p1": ["A", "B"], "p2": ["C"], "p3": ["A"]}
    sim = spatial_similarity(SAMPLES, normalize=False, pools=pools)
    assert list(sim.index) == ["p1", "p2", "p3"]
    assert_sim(sim, 0, 1 - np.sqrt(6) / np.sqrt(38), 1 - 4 / np.sqrt(38))
