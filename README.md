# The Laplacian in geostatistics

Interactive explainer: <https://gstvschlz.com/laplacian/>

One animation of what the Laplacian is and what it does in geostatistical
estimation and simulation: [`laplacian.mp4`](laplacian.mp4), 57 seconds.

1. **The operator.** The Laplacian of a cell is its value minus the average
   of its neighbours: positive on a high, negative in a low, zero on a flat
   or a straight slope.
2. **Estimation.** Hold the data and make the Laplacian zero everywhere
   else. Replacing each cell by the average of its neighbours, repeatedly,
   gets there: straight lines between the data, which is kriging with a
   linear variogram.
3. **κ.** Divide by a little more than the number of neighbours and the
   estimate sags to the mean away from the data. `κ` sets the range; this is
   simple kriging with an exponential variogram.
4. **On a map.** The same rule with four neighbours gives the kriging map.
5. **Simulation.** Feed the same operator with white noise, solving
   `(κ² I + L) z = noise`. The result is a field with a Matérn variogram.
6. **Both.** Hold the data and feed noise: conditional simulation. Without
   the noise it is kriging.

The model behind all of it is a Gaussian field with precision (inverse
covariance) `Q = τ(κ² I + L)^α`, with `L` the Laplacian: `κ` is the range,
`α` the smoothness at the origin, `τ` the sill (Lindgren, Rue and Lindström,
2011).

## Run it

```sh
npm run serve      # the interactive explainer at http://localhost:8000
npm test           # runs the Node.js unit checks in test/
mise run animate   # rewrites laplacian.mp4, about two minutes; needs ffmpeg
mise run verify    # lint, then the theory checks
```

## Layout

- `index.html`, `style.css`, `web/` — the six-step explainer, no build step;
  published by GitHub Pages from `main`. Theme and fonts follow gstvschlz.com.
- `test/` — Node.js unit checks for `web/laplacian.js` and `web/walker.js`.
- `animate.py` — the six scenes.
- `src/lapgeo/` — the Laplacian, the precision, draws, conditioning and the
  cellular rule the scenes use.
- `tests/test_theory.py` — checks of what the animation states, with seeded
  draws.
