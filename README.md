# The Laplacian in geostatistics

What the Laplacian is, and what it is doing inside a kriging system or a
simulation, worked out on a grid with code that checks every statement.

## The short answer

The Laplacian of a field at a cell is **the value of the cell minus the
average of its neighbours** (times the number of neighbours). It is zero where
the field is flat or a ramp, positive on a high, negative in a low.

In geostatistics it is not used to measure fields. It is used to *define*
them, in place of a variogram, and it shows up in three roles that are the
same statement read three ways:

1. **As kriging weights.** Write the model as a precision matrix (the inverse
   of the covariance matrix) `Q = τ(κ² I + L)^α`, with `L` the Laplacian. Then
   the simple kriging of a cell from *all* the other cells puts weight only on
   the stencil of `L`: four neighbours get `1/(4 + κ²)` each and every other
   cell gets exactly zero. The Laplacian is a screen effect that is exact.
2. **As a variogram model.** That precision is a Matérn covariance of
   smoothness `ν = α − d/2` in `d` dimensions, with `κ` setting the range and
   `τ` the sill. A tridiagonal precision in 1-D is the exponential model. In
   3-D the exponential model is `α = 2`.
3. **As an equation for the estimate.** The kriging mean solves
   `(κ² I + L)^α m = 0` between the data, with the data as boundary values.
   With `κ = 0` and `α = 1` that is Laplace's equation: every estimated cell
   is the average of its neighbours. With `α = 2` it is a spline.

Simulation follows from the same object: filter white noise through the
inverse of `(κ² I + L)^(α/2)`, which for `α = 2` is one sparse solve of
`(κ² I + L) z = w`. This is the stochastic partial differential equation
(SPDE) approach of Lindgren, Rue and Lindström (2011).

The parameters, in variogram terms:

| Symbol | Role | Variogram equivalent |
| --- | --- | --- |
| `L` | cell minus neighbours | the shape of the model |
| `κ` | pull of each cell towards the mean | range ≈ `3/κ` (ν = ½), `4/κ` (ν = 1), at 95 % of the sill |
| `α` | how many times the operator is applied | behaviour at the origin, `ν = α − d/2` |
| `τ` | scale of the precision | sill `∝ 1/(τ κ^2ν)` |
| link weights | strength of each neighbour | anisotropy, barriers |

## Chapters

Each chapter is a short text, one figure, a script that draws it and tests
that hold it to account.

1. [Neighbours](chapters/01-neighbours/) — the operator itself.
2. [The screen effect](chapters/02-screen-effect/) — precision against covariance, and why the weights are the stencil.
3. [Which variogram](chapters/03-which-variogram/) — Matérn, `ν = α − d/2`, range and sill.
4. [Simulation](chapters/04-simulation/) — three ways to draw, and conditioning.
5. [Estimation](chapters/05-estimation/) — kriging as a differential equation, splines, overshoot.
6. [No sill](chapters/06-no-sill/) — `κ → 0`, the intrinsic hypothesis and ordinary kriging.
7. [Anisotropy and barriers](chapters/07-anisotropy-and-barriers/) — changing the links, and what the approach cannot do.

## Run it

```sh
mise run test      # the theory checks, about half a minute
mise run figures   # redraws every chapters/*/*.png
mise run verify    # lint, then test
mise run preview   # one page with every chapter, on http://localhost:8791/_site/
```

Without mise: `uv sync`, then `uv run pytest` and
`uv run python chapters/01-neighbours/run.py`.

## Layout

- `src/lapgeo/__init__.py` — the Laplacian, the precision, exact lattice
  covariances, draws, conditioning, the cellular rule and the Matérn model.
- `src/lapgeo/style.py` — the figure style.
- `chapters/NN-name/` — `README.md`, `run.py` and the figure it writes.
- `tests/test_theory.py` — one check per statement, grouped by chapter, with
  seeded draws.

## Conventions

- Cells have unit size; lags are in cells. `L` is the positive semi-definite
  graph Laplacian, the discrete `−Δ`, so it is positive on a high.
- Fields have mean zero. Data sit on cells and are honoured exactly; there is
  no nugget.
- Most grids are tori (they wrap around), which makes the Fourier modes exact
  eigenvectors and the checks exact.

## References

- Whittle, P. (1954). On stationary processes in the plane. *Biometrika* 41.
- Matheron, G. (1973). The intrinsic random functions and their applications. *Advances in Applied Probability* 5.
- Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. *JRSS B* 36.
- Matheron, G. (1981). Splines and kriging: their formal equivalence. In *Down-to-Earth Statistics*, Syracuse University.
- Rue, H. and Held, L. (2005). *Gaussian Markov Random Fields*. Chapman & Hall.
- Lindgren, F., Rue, H. and Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *JRSS B* 73.
