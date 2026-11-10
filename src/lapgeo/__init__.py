"""The lattice Laplacian as a geostatistical model.

Cells are unit squares (or segments, or cubes). ``L`` is the positive
semi-definite graph Laplacian, the discrete version of ``-Δ``: applied to a
field it returns, at every cell, the number of neighbours times the value
minus the sum of the neighbours. The model is a Gaussian field with precision
(inverse covariance) ``Q = τ (κ² I + L)^α``.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.special import gamma, kv

__all__ = [
    "conditional_mean",
    "conditional_simulation",
    "covariance",
    "eigenvalues",
    "grid_edges",
    "laplacian",
    "laplacian_from_edges",
    "matern",
    "matern_variance",
    "precision",
    "relax",
    "simple_kriging",
    "simulate_spde",
    "simulate_spectral",
]


def _weights(shape, weights):
    return (1.0,) * len(shape) if weights is None else tuple(float(w) for w in weights)


def _second_difference(m, periodic):
    off = -np.ones(m - 1)
    d = sp.diags_array([off, off], offsets=[-1, 1], shape=(m, m), format="lil")
    if periodic:
        d[0, m - 1] += -1.0
        d[m - 1, 0] += -1.0
    d = d.tocsr()
    return sp.diags_array(-d.sum(axis=1)) + d


def laplacian(shape, periodic=True, weights=None):
    """Graph Laplacian of a regular grid.

    Parameters
    ----------
    shape : tuple of int
        Cells along each axis; the field is flattened in C order.
    periodic : bool
        Wrap the grid into a torus. Otherwise border cells have fewer
        neighbours, which is a no-flux (reflecting) boundary.
    weights : sequence of float, optional
        One weight per axis on the links along it. Unequal weights give
        geometric anisotropy; the range along an axis grows with the square
        root of its weight.

    Returns
    -------
    scipy.sparse.csr_array
        Symmetric, positive semi-definite, rows summing to zero.
    """
    shape = tuple(shape)
    total = sp.csr_array((int(np.prod(shape)),) * 2)
    for axis, (m, w) in enumerate(zip(shape, _weights(shape, weights))):
        term = sp.eye_array(1, format="csr")
        for other, k in enumerate(shape):
            term = sp.kron(term, _second_difference(m, periodic) if other == axis else sp.eye_array(k), "csr")
        total = total + w * term
    return total.tocsr()


def grid_edges(shape):
    """Pairs of flat indices of neighbouring cells of a non-periodic grid."""
    index = np.arange(int(np.prod(shape))).reshape(shape)
    pairs = []
    for axis in range(len(shape)):
        a = np.take(index, np.arange(shape[axis] - 1), axis=axis).ravel()
        b = np.take(index, np.arange(1, shape[axis]), axis=axis).ravel()
        pairs.append(np.column_stack([a, b]))
    return np.vstack(pairs)


def laplacian_from_edges(n, edges, weights=None):
    """Graph Laplacian ``D - A`` of any set of links between ``n`` cells."""
    edges = np.asarray(edges)
    w = np.ones(len(edges)) if weights is None else np.asarray(weights, dtype=float)
    i, j = edges[:, 0], edges[:, 1]
    a = sp.coo_array((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(n, n)).tocsr()
    return (sp.diags_array(a.sum(axis=1)) - a).tocsr()


def precision(shape=None, kappa=0.0, alpha=1, tau=1.0, periodic=True, weights=None, lap=None):
    """Precision matrix ``τ (κ² I + L)^α`` for integer ``α``."""
    lap = laplacian(shape, periodic, weights) if lap is None else lap
    a = (kappa**2 * sp.eye_array(lap.shape[0]) + lap).tocsr()
    q = a
    for _ in range(int(alpha) - 1):
        q = q @ a
    return (tau * q).tocsr()


def eigenvalues(shape, weights=None):
    """Eigenvalues of the torus Laplacian, laid out like ``numpy.fft.fftn``.

    The eigenvectors are the Fourier modes, and the eigenvalue of the mode of
    angular frequency ``ω`` per cell is ``Σ w (2 - 2 cos ω)``, close to
    ``|ω|²`` for long waves.
    """
    shape = tuple(shape)
    lam = np.zeros(shape)
    for axis, (m, w) in enumerate(zip(shape, _weights(shape, weights))):
        one = w * (2.0 - 2.0 * np.cos(2.0 * np.pi * np.fft.fftfreq(m)))
        lam = lam + one.reshape([-1 if a == axis else 1 for a in range(len(shape))])
    return lam


def _spectrum(shape, kappa, alpha, tau, weights):
    base = kappa**2 + eigenvalues(shape, weights)
    with np.errstate(divide="ignore"):
        return np.where(base > 0, 1.0 / (tau * base**alpha), 0.0)


def covariance(shape, kappa, alpha=1, tau=1.0, weights=None):
    """Exact covariance of the torus model between cell 0 and every cell.

    ``α`` may be any positive real. With ``κ = 0`` the constant mode is
    dropped (a pseudo-inverse): the values are then only meaningful through
    differences, ``γ(h) = C(0) - C(h)``.
    """
    return np.fft.ifftn(_spectrum(shape, kappa, alpha, tau, weights)).real


def simulate_spectral(shape, kappa, alpha=1, tau=1.0, rng=None, weights=None):
    """Unconditional draw: white noise filtered by ``Q^(-1/2)`` in Fourier space."""
    rng = np.random.default_rng(rng)
    root = np.sqrt(_spectrum(shape, kappa, alpha, tau, weights))
    return np.fft.ifftn(root * np.fft.fftn(rng.standard_normal(shape))).real


def simulate_spde(shape, kappa, tau=1.0, rng=None, periodic=True, weights=None, lap=None):
    """Unconditional draw for ``α = 2`` by solving ``(κ² I + L) z = w / √τ``.

    This is the discretised stochastic partial differential equation
    ``(κ² - Δ) Z = W``: one sparse solve with white noise on the right.
    """
    rng = np.random.default_rng(rng)
    a = precision(shape, kappa, 1, 1.0, periodic, weights, lap)
    n = a.shape[0]
    z = splu(a.tocsc()).solve(rng.standard_normal(n) / np.sqrt(tau))
    return z if shape is None else z.reshape(shape)


def conditional_mean(q, data_index, data_values):
    """Mean of a zero-mean field with precision ``q`` given exact data.

    Solves ``Q_uu m_u = -Q_ud z_d`` over the unknown cells: the rows of the
    precision at the unknown cells, applied to the full field, must be zero.
    This is simple kriging with a global neighbourhood.
    """
    n = q.shape[0]
    data_index = np.asarray(data_index)
    known = np.zeros(n, dtype=bool)
    known[data_index] = True
    u = np.flatnonzero(~known)
    q = q.tocsr()
    rows = q[u]
    out = np.empty(n)
    out[data_index] = data_values
    out[u] = splu(rows[:, u].tocsc()).solve(-(rows[:, data_index] @ np.asarray(data_values, dtype=float)))
    return out


def conditional_simulation(q, data_index, data_values, unconditional):
    """Condition an unconditional draw of the same model by kriging its error."""
    z = np.asarray(unconditional, dtype=float).ravel()
    data_index = np.asarray(data_index)
    return z + conditional_mean(q, data_index, np.asarray(data_values) - z[data_index])


def simple_kriging(cov, data_index, data_values):
    """Dense simple kriging of every cell, mean zero, global neighbourhood.

    Returns the estimates and the weights (one row per cell).
    """
    data_index = np.asarray(data_index)
    weights = np.linalg.solve(cov[np.ix_(data_index, data_index)], cov[data_index]).T
    return weights @ np.asarray(data_values, dtype=float), weights


def relax(shape, kappa, data_index, data_values, sweeps, tau=1.0, rng=None, start=None, source=None):
    """Cellular rule on a 2-D torus with ``α = 1``.

    Each sweep replaces every free cell by the sum of its four neighbours
    divided by ``κ² + 4``, one colour of the chequerboard at a time. Without
    ``rng`` the rule converges to the kriging mean (it is Gauss-Seidel on
    ``Q_uu m = -Q_ud z_d``). With ``rng`` it adds noise of variance
    ``1 / (τ (κ² + 4))`` and becomes a Gibbs sampler of the conditional field.
    With ``source`` (one fixed number per cell, added to the sum of the
    neighbours) it converges to the solution of ``(κ² I + L) z = source``.
    """
    ny, nx = shape
    if ny % 2 or nx % 2:
        raise ValueError("the chequerboard needs even sides on a torus")
    rng = None if rng is None else np.random.default_rng(rng)
    z = np.zeros(shape) if start is None else np.array(start, dtype=float)
    data_index = np.asarray(data_index, dtype=int)
    free = np.ones(shape, dtype=bool)
    free.flat[data_index] = False
    z.flat[data_index] = data_values
    iy, ix = np.indices(shape)
    diag = kappa**2 + 4.0
    for _ in range(sweeps):
        for colour in (0, 1):
            mask = free & ((iy + ix) % 2 == colour)
            nb = np.roll(z, 1, 0) + np.roll(z, -1, 0) + np.roll(z, 1, 1) + np.roll(z, -1, 1)
            new = nb / diag if source is None else (nb + source) / diag
            if rng is not None:
                new = new + rng.standard_normal(shape) / np.sqrt(tau * diag)
            z[mask] = new[mask]
    return z


def matern_variance(kappa, alpha, d, tau=1.0):
    """Sill of the continuous field with precision ``τ (κ² - Δ)^α`` in ``d`` dimensions."""
    nu = alpha - d / 2.0
    return gamma(nu) / (gamma(alpha) * (4.0 * np.pi) ** (d / 2.0) * kappa ** (2.0 * nu) * tau)


def matern(h, kappa, nu, sill=1.0):
    """Matérn covariance ``sill · 2^(1-ν)/Γ(ν) · (κh)^ν K_ν(κh)``."""
    x = kappa * np.abs(np.asarray(h, dtype=float))
    out = np.full(x.shape, float(sill))
    pos = x > 0
    out[pos] = sill * 2.0 ** (1.0 - nu) / gamma(nu) * x[pos] ** nu * kv(nu, x[pos])
    return out
