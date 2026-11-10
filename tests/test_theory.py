"""The statements the animation makes, and the theory behind them, as checks that can be rerun."""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

import lapgeo as lg


def dense_covariance(shape, kappa, alpha):
    return np.linalg.inv(lg.precision(shape, kappa, alpha).toarray())


# 01 — neighbours


def test_laplacian_is_cell_minus_neighbours():
    rng = np.random.default_rng(0)
    z = rng.standard_normal((12, 10))
    neighbours = np.roll(z, 1, 0) + np.roll(z, -1, 0) + np.roll(z, 1, 1) + np.roll(z, -1, 1)
    assert np.allclose(lg.laplacian(z.shape) @ z.ravel(), (4 * z - neighbours).ravel())


def test_laplacian_ignores_constants_and_ramps_and_sees_curvature():
    shape = (9, 11)
    y, x = np.indices(shape)
    lap = lg.laplacian(shape, periodic=False)
    inner = (slice(1, -1), slice(1, -1))
    assert np.allclose(lap @ np.ones(99), 0)
    assert np.allclose((lap @ (3 * x - 2 * y + 1.0).ravel()).reshape(shape)[inner], 0)
    assert np.allclose((lap @ (x**2.0).ravel()).reshape(shape)[inner], -2)


def test_waves_are_eigenvectors():
    shape = (16, 20)
    y, x = np.indices(shape)
    ky, kx = 3, 5
    wave = np.cos(2 * np.pi * (ky * y / shape[0] + kx * x / shape[1]))
    factor = 4 - 2 * np.cos(2 * np.pi * ky / shape[0]) - 2 * np.cos(2 * np.pi * kx / shape[1])
    assert np.allclose(lg.laplacian(shape) @ wave.ravel(), factor * wave.ravel())
    assert np.isclose(lg.eigenvalues(shape)[ky, kx], factor)


# 02 — screen effect


def test_one_dimensional_precision_is_the_exponential_covariance():
    kappa, n = 0.3, 400
    c = lg.covariance((n,), kappa, 1)
    rho = 1 + kappa**2 / 2 - np.sqrt(kappa**2 + kappa**4 / 4)
    h = np.arange(60)
    assert np.allclose(c[h], c[0] * rho**h, rtol=1e-10)
    assert np.isclose(c[0], 1 / (kappa * np.sqrt(4 + kappa**2)))
    assert np.allclose(c, dense_covariance((n,), kappa, 1)[0])


def test_kriging_weights_are_the_stencil():
    shape, kappa = (14, 14), 0.5
    cov = dense_covariance(shape, kappa, 1)
    others = np.arange(1, cov.shape[0])
    _, w = lg.simple_kriging(cov, others, np.zeros(others.size))
    weights = np.r_[0.0, w[0]].reshape(shape)
    expected = np.zeros(shape)
    expected[0, 1] = expected[1, 0] = expected[0, -1] = expected[-1, 0] = 1 / (4 + kappa**2)
    assert np.allclose(weights, expected, atol=1e-12)
    variance = cov[0, 0] - w[0] @ cov[others, 0]
    assert np.isclose(variance, 1 / (4 + kappa**2))


# 03 — which variogram


@pytest.mark.parametrize(
    ("shape", "alpha", "tolerance"),
    [((4096,), 1, 0.01), ((512, 512), 2, 0.03), ((96, 96, 96), 2, 0.07)],
)
def test_lattice_covariance_is_matern(shape, alpha, tolerance):
    kappa, d = 0.2, len(shape)
    h = np.arange(0, 17)
    lattice = lg.covariance(shape, kappa, alpha)[(0,) * (d - 1)][h]
    model = lg.matern(h, kappa, alpha - d / 2, lg.matern_variance(kappa, alpha, d))
    assert np.allclose(lattice, model, rtol=tolerance)


def test_two_dimensional_alpha_one_has_no_sill_of_its_own():
    sills = [lg.covariance((1024, 1024), kappa, 1)[0, 0] for kappa in (0.2, 0.1, 0.05)]
    assert np.allclose(np.diff(sills), np.log(2) / (2 * np.pi), rtol=0.02)


# 04 — simulation


def test_conditioning_with_the_precision_is_simple_kriging():
    rng = np.random.default_rng(1)
    shape = (18, 18)
    data = rng.choice(324, 12, replace=False)
    values = rng.standard_normal(12)
    for alpha in (1, 2):
        q = lg.precision(shape, 0.3, alpha)
        kriged, _ = lg.simple_kriging(np.linalg.inv(q.toarray()), data, values)
        assert np.allclose(lg.conditional_mean(q, data, values), kriged, atol=1e-10)


def test_draws_reproduce_the_variogram_and_the_data():
    shape, kappa = (128, 128), 0.25
    tau = lg.matern_variance(kappa, 2, 2)
    rng = np.random.default_rng(4)
    exact = lg.covariance(shape, kappa, 2, tau)
    h = np.array([1, 2, 4, 8, 16])
    exact_gamma = exact[0, 0] - exact[0, h]
    for draw in (lg.simulate_spde, lambda s, k, t, r: lg.simulate_spectral(s, k, 2, t, r)):
        fields = np.stack([draw(shape, kappa, tau, rng) for _ in range(40)])
        gamma = [0.5 * np.mean((fields - np.roll(fields, lag, axis=2)) ** 2) for lag in h]
        assert np.allclose(gamma, exact_gamma, rtol=0.06)

    q = lg.precision(shape, kappa, 2, tau)
    data = rng.choice(shape[0] * shape[1], 25, replace=False)
    values = rng.standard_normal(25)
    z = lg.conditional_simulation(q, data, values, lg.simulate_spde(shape, kappa, tau, rng))
    assert np.allclose(z[data], values)


def test_gibbs_rule_samples_the_model():
    shape, kappa = (32, 32), 0.6
    rng = np.random.default_rng(8)
    z = np.zeros(shape)
    total = np.zeros(3)
    kept = 0
    for sweep in range(3000):
        z = lg.relax(shape, kappa, [], [], 1, rng=rng, start=z)
        if sweep >= 200:
            total += [np.mean(z * np.roll(z, lag, axis=1)) for lag in (0, 1, 2)]
            kept += 1
    assert np.allclose(total / kept, lg.covariance(shape, kappa, 1)[0, :3], rtol=0.05)


# 05 — estimation


def test_kriging_mean_minimises_the_energy():
    rng = np.random.default_rng(2)
    shape = (16, 16)
    data = rng.choice(256, 10, replace=False)
    values = rng.standard_normal(10)
    for alpha in (1, 2):
        q = lg.precision(shape, 0.2, alpha)
        m = lg.conditional_mean(q, data, values)
        for _ in range(20):
            other = m + 0.1 * rng.standard_normal(256)
            other[data] = values
            assert other @ (q @ other) > m @ (q @ m)


def test_relaxation_converges_to_the_kriging_mean():
    rng = np.random.default_rng(3)
    shape = (24, 24)
    data = rng.choice(576, 15, replace=False)
    values = rng.standard_normal(15)
    q = lg.precision(shape, 0.3, 1)
    relaxed = lg.relax(shape, 0.3, data, values, 400)
    assert np.allclose(relaxed.ravel(), lg.conditional_mean(q, data, values), atol=1e-10)


def test_rule_fed_with_noise_solves_the_spde():
    shape, kappa = (32, 32), 0.3
    noise = np.random.default_rng(9).standard_normal(shape)
    relaxed = lg.relax(shape, kappa, [], [], 600, source=noise)
    assert np.allclose(lg.precision(shape, kappa, 1) @ relaxed.ravel(), noise.ravel(), atol=1e-9)


def test_mixing_two_draws_by_cos_and_sin_keeps_the_data():
    rng = np.random.default_rng(10)
    shape, kappa = (24, 24), 0.3
    q = lg.precision(shape, kappa, 2)
    data = rng.choice(576, 10, replace=False)
    values = rng.standard_normal(10)
    mean = lg.conditional_mean(q, data, values)
    draws = [lg.simulate_spectral(shape, kappa, 2, rng=rng).ravel() for _ in range(2)]
    residuals = [z - lg.conditional_mean(q, data, z[data]) for z in draws]
    mixed = mean + np.cos(0.7) * residuals[0] + np.sin(0.7) * residuals[1]
    assert np.allclose(mixed[data], values)


def test_alpha_one_never_overshoots_and_alpha_two_does():
    rng = np.random.default_rng(1)
    shape = (20, 20)
    data = rng.choice(400, 12, replace=False)
    values = rng.standard_normal(12)
    q = lg.precision(shape, 0.1, 1)
    m = lg.conditional_mean(q, data, values)
    _, w = lg.simple_kriging(np.linalg.inv(q.toarray()), data, values)
    assert w.min() > -1e-12
    assert values.min() - 1e-12 <= m.min() and m.max() <= values.max() + 1e-12

    pair = np.ravel_multi_index(([10, 10], [8, 11]), shape)
    q = lg.precision(shape, 0.1, 2)
    m = lg.conditional_mean(q, pair, [1.0, 2.0])
    _, w = lg.simple_kriging(np.linalg.inv(q.toarray()), pair, [1.0, 2.0])
    assert w.min() < -0.05
    assert m.max() > 2.05


def test_one_dimension_gives_linear_interpolation_then_the_cubic_spline():
    n = 240
    data = np.array([20, 55, 75, 130, 150, 205])
    values = np.array([0.2, 1.0, 1.6, -0.4, -0.9, 0.5])
    x = np.arange(n)
    ring_x = np.r_[data, data[0] + n]
    ring_v = np.r_[values, values[0]]
    wrapped = np.where(x < data[0], x + n, x)

    linear = lg.conditional_mean(lg.precision((n,), 0.0, 1), data, values)
    assert np.allclose(linear, np.interp(wrapped, ring_x, ring_v))

    cubic = lg.conditional_mean(lg.precision((n,), 0.0, 2), data, values)
    spline = CubicSpline(ring_x, ring_v, bc_type="periodic")(wrapped)
    assert np.abs(cubic - spline).max() < 0.01


# 06 — no sill


def test_kappa_squared_is_the_weight_of_the_mean():
    shape = (12, 12)
    for kappa in (1.0, 0.3, 0.05):
        cov = dense_covariance(shape, kappa, 1)
        others = np.arange(1, cov.shape[0])
        _, w = lg.simple_kriging(cov, others, np.zeros(others.size))
        assert np.isclose(w[0].sum(), 4 / (4 + kappa**2))


def test_limit_variograms():
    n = 200
    c = lg.covariance((n,), 0.0, 1)
    h = np.arange(n)
    assert np.allclose(c[0] - c, h * (n - h) / (2 * n))

    lags = np.array([4, 8, 16, 32, 64])
    c = lg.covariance((1024, 1024), 0.0, 1)
    assert np.allclose(np.diff(c[0, 0] - c[0, lags]), np.log(2) / (2 * np.pi), rtol=0.01)

    c = lg.covariance((1024, 1024), 0.0, 2)
    scaled = (c[0, 0] - c[0, lags[1:]]) / lags[1:] ** 2
    assert np.allclose(np.diff(scaled), -np.log(2) / (8 * np.pi), rtol=0.05)


def test_intrinsic_estimate_does_not_see_a_shift_of_the_data():
    rng = np.random.default_rng(6)
    shape = (16, 16)
    data = rng.choice(256, 9, replace=False)
    values = rng.standard_normal(9)
    q = lg.precision(shape, 0.0, 1)
    shifted = lg.conditional_mean(q, data, values + 100.0)
    assert np.allclose(shifted, lg.conditional_mean(q, data, values) + 100.0)


# 07 — anisotropy and barriers


def test_link_weights_stretch_the_range_by_their_square_root():
    c = lg.covariance((512, 512), 0.2, 2, weights=(1.0, 4.0))
    h = np.array([2, 4, 8])
    assert np.allclose(c[0, 2 * h], c[h, 0], rtol=0.01)


def test_cutting_links_cuts_covariance():
    shape = (30, 30)
    n = 900
    edges = lg.grid_edges(shape)
    _, ix = np.unravel_index(edges, shape)
    crosses = (ix.min(axis=1) == 14) & (ix.max(axis=1) == 15)
    whole = lg.laplacian_from_edges(n, edges)
    assert abs(whole - lg.laplacian(shape, periodic=False)).max() < 1e-12
    cut = lg.laplacian_from_edges(n, edges[~crosses])
    a, b = np.ravel_multi_index(([15, 15], [13, 16]), shape)
    with_link = np.linalg.inv(lg.precision(kappa=0.2, alpha=2, lap=whole).toarray())
    without = np.linalg.inv(lg.precision(kappa=0.2, alpha=2, lap=cut).toarray())
    assert with_link[a, b] > 0.3 * with_link[a, a]
    assert abs(without[a, b]) < 1e-12
