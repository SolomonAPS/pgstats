"""Neutral-coalescent calibration tests for standardized neutrality statistics.

Under the standard neutral model with no recombination (the model the variance
formulas are derived under), Tajima's D, Fu & Li's D/D*/F/F*, Zeng's E and the
normalized Fay & Wu's H must all come out with mean ~0 and sd ~1, and produce
no NaNs for polymorphic windows.

This is the test that catches coefficient errors the unit tests cannot: the
b_n-vs-b2 notation collision (negative Var(D*) -> NaNs above a sample-size-
dependent S threshold), Tajima's missing e1/e2 step (~a1-fold deflation), the
theta_H factor-2 omission (+theta/2 bias in H), and any numerator/variance
convention mismatch. The (n=74, S~400) cell is the one that catches the
negative-variance regime seen in real data.
"""
import numpy as np
import pytest
import xarray as xr

from pgstats.stats.sfs_statistics import (
    tajima_d, fu_li_d, fu_li_f, zeng_e, fay_wu_h, calculate_a1,
)

REPS = 200
CONFIGS = [(20, 30), (74, 100), (74, 400)]
SEED = 20260812

STATS = [
    ("tajima_d",     lambda ds: tajima_d(ds)["tajima_d"].values),
    ("fu_li_d_star", lambda ds: fu_li_d(ds, folded=True)["fu_li_d_star"].values),
    ("fu_li_d",      lambda ds: fu_li_d(ds, folded=False)["fu_li_d"].values),
    ("fu_li_f_star", lambda ds: fu_li_f(ds, folded=True)["fu_li_f_star"].values),
    ("fu_li_f",      lambda ds: fu_li_f(ds, folded=False)["fu_li_f"].values),
    ("zeng_e",       lambda ds: zeng_e(ds)["zeng_e"].values),
    ("fay_wu_h",     lambda ds: fay_wu_h(ds)["fay_wu_h"].values),
]


def _sim_matrix(n, theta, rng):
    """One neutral Kingman coalescent replicate -> (S, n) 0/1 variant matrix."""
    active = [(1 << i, 0.0) for i in range(n)]
    branches = []
    k = n
    while k > 1:
        rate = k * (k - 1) / 2.0
        dt = rng.exponential(1.0 / rate)
        active = [(m, l + dt) for (m, l) in active]
        i, j = sorted(rng.choice(k, size=2, replace=False))
        branches.append(active[i])
        branches.append(active[j])
        merged = (active[i][0] | active[j][0], 0.0)
        active = [active[x] for x in range(k) if x != i and x != j] + [merged]
        k -= 1
    masks = []
    for mask, length in branches:
        for _ in range(rng.poisson(theta / 2.0 * length)):
            masks.append(mask)
    mat = np.zeros((len(masks), n), dtype=np.int8)
    for r, mask in enumerate(masks):
        for b in range(n):
            if (mask >> b) & 1:
                mat[r, b] = 1
    return mat


def _build_dataset(n, theta, reps, rng):
    mats = [_sim_matrix(n, theta, rng) for _ in range(reps)]
    sizes = [m.shape[0] for m in mats]
    stops = np.cumsum(sizes).astype(np.int64)
    starts = np.concatenate([[0], stops[:-1]]).astype(np.int64)
    geno = np.concatenate(mats)[:, :, None]
    return xr.Dataset(
        {
            "call_genotype": (("variants", "samples", "ploidy"), geno),
            "window_start_idx": ("windows", starts),
            "window_stop_idx": ("windows", stops),
        },
        coords={"windows": np.arange(reps)},
    )


@pytest.fixture(scope="module", params=CONFIGS, ids=lambda c: f"n{c[0]}_S{c[1]}")
def neutral_dataset(request):
    n, s_target = request.param
    rng = np.random.default_rng(SEED + n + s_target)
    theta = s_target / calculate_a1(n)
    return _build_dataset(n, theta, REPS, rng)


@pytest.mark.parametrize("stat_name,stat_fn", STATS, ids=[s[0] for s in STATS])
def test_neutral_calibration(neutral_dataset, stat_name, stat_fn):
    vals = stat_fn(neutral_dataset)
    finite = np.isfinite(vals)
    frac_nan = 1 - finite.mean()
    v = vals[finite]
    assert frac_nan < 0.02, f"{stat_name}: {100 * frac_nan:.1f}% NaN windows"
    assert abs(v.mean()) < 0.35, f"{stat_name}: mean {v.mean():.2f} (expected ~0)"
    assert 0.6 < v.std() < 1.3, f"{stat_name}: sd {v.std():.2f} (expected ~1)"
    frac_extreme = np.mean(np.abs(v) > 4)
    assert frac_extreme < 0.025, (
        f"{stat_name}: {100 * frac_extreme:.1f}% of values beyond |4| "
        "(exploding normalization)")
