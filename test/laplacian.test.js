import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import {
  laplacian1D,
  discrepancy1D,
  laplacian2D,
  relax1D,
  solve1D,
  relax2D,
  conditionalSimulate2D,
  mulberry32
} from "../web/laplacian.js";
import { parseCsv, makeCase } from "../web/walker.js";

test("Laplacian of a linear field is 0 on interior cells (1-D and 2-D)", () => {
  const n = 25;
  const line = Float64Array.from({ length: n }, (_, i) => 1.7 * i - 3.2);
  const lap1 = laplacian1D(line, { periodic: false });
  const disc1 = discrepancy1D(line, { periodic: false });
  for (let i = 1; i < n - 1; i++) {
    assert.ok(Math.abs(lap1[i]) < 1e-12, `1D interior cell ${i} has non-zero Laplacian ${lap1[i]}`);
    assert.ok(Math.abs(disc1[i]) < 1e-12, `1D interior cell ${i} has non-zero discrepancy ${disc1[i]}`);
  }

  // Peak is positive, pit is negative
  const bump = new Float64Array(9);
  bump[2] = 1.0;
  bump[6] = -1.0;
  const bumpLap = laplacian1D(bump);
  assert.ok(bumpLap[2] > 0, "Peak should have positive Laplacian");
  assert.ok(bumpLap[6] < 0, "Pit should have negative Laplacian");

  // 2-D linear ramp
  const w = 12, h = 10;
  const plane = new Float64Array(w * h);
  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      plane[y * w + x] = 2.5 * x - 1.1 * y + 4.0;
    }
  }
  const lap2 = laplacian2D(plane, w, h, { periodic: false });
  for (let y = 1; y < h - 1; y++) {
    for (let x = 1; x < w - 1; x++) {
      assert.ok(Math.abs(lap2[y * w + x]) < 1e-12, `2D interior cell (${x},${y}) has non-zero Laplacian`);
    }
  }
});

test("1-D relaxation converges to linear interpolation between data points", () => {
  const n = 33;
  const dataIndices = [4, 14, 26];
  const dataValues = [0.5, 2.0, -1.0];
  const z = new Float64Array(n);
  relax1D(z, { kappa: 0, dataIndices, dataValues, sweeps: 600, periodic: false });

  // Outside first datum: flat at dataValues[0] (reflecting boundary)
  for (let i = 0; i <= dataIndices[0]; i++) {
    assert.ok(Math.abs(z[i] - dataValues[0]) < 1e-5, `Cell ${i} should match left boundary datum`);
  }
  // Between data points: linear segments
  for (let seg = 0; seg < dataIndices.length - 1; seg++) {
    const i0 = dataIndices[seg], i1 = dataIndices[seg + 1];
    const v0 = dataValues[seg], v1 = dataValues[seg + 1];
    for (let i = i0; i <= i1; i++) {
      const t = (i - i0) / (i1 - i0);
      const expected = v0 + t * (v1 - v0);
      assert.ok(Math.abs(z[i] - expected) < 1e-4, `Cell ${i} expected ${expected}, got ${z[i]}`);
    }
  }
  // Outside last datum: flat at dataValues[last]
  for (let i = dataIndices[2]; i < n; i++) {
    assert.ok(Math.abs(z[i] - dataValues[2]) < 1e-5, `Cell ${i} should match right boundary datum`);
  }

  // Also matches exact tridiagonal solver solve1D
  const exact = solve1D(n, { kappa: 0, dataIndices, dataValues });
  for (let i = 0; i < n; i++) {
    assert.ok(Math.abs(z[i] - exact[i]) < 1e-4);
  }
});

test("kappa > 0 sags the estimate toward the mean away from data", () => {
  const n = 41;
  const dataIndices = [5, 35];
  const dataValues = [1.5, 1.5];

  const linear = solve1D(n, { kappa: 0, dataIndices, dataValues });
  const sagMild = solve1D(n, { kappa: 0.15, dataIndices, dataValues });
  const sagStrong = solve1D(n, { kappa: 0.4, dataIndices, dataValues });

  const mid = 20;
  assert.ok(Math.abs(linear[mid] - 1.5) < 1e-9, "kappa = 0 stays flat at 1.5 between equal data");
  assert.ok(sagMild[mid] < linear[mid] - 0.3, "kappa = 0.15 sags toward 0");
  assert.ok(sagStrong[mid] < sagMild[mid] - 0.15, "larger kappa sags closer to 0");
  assert.ok(sagStrong[mid] > 0 && sagStrong[mid] < 0.05, "stays positive and close to 0 mean");

  // Relaxation with kappa > 0 converges to solve1D
  const relaxed = new Float64Array(n);
  relax1D(relaxed, { kappa: 0.25, dataIndices, dataValues, sweeps: 300 });
  const exact = solve1D(n, { kappa: 0.25, dataIndices, dataValues });
  for (let i = 0; i < n; i++) {
    assert.ok(Math.abs(relaxed[i] - exact[i]) < 1e-6);
  }
});

test("conditional simulation reproduces held data and collapses to estimate when noise is off", () => {
  const w = 24, h = 24;
  const dataIndices = [2 * w + 3, 10 * w + 12, 18 * w + 7, 20 * w + 20];
  const dataValues = [1.8, -1.2, 0.9, -0.5];
  const rng = mulberry32(42);
  const noise = rng.normalArray(w * h);

  const sim = conditionalSimulate2D(w, h, {
    kappa: 0.2,
    dataIndices,
    dataValues,
    source: noise,
    noiseScale: 1.0,
    sweeps: 300
  });

  const est = conditionalSimulate2D(w, h, {
    kappa: 0.2,
    dataIndices,
    dataValues,
    source: noise,
    noiseScale: 0.0,
    sweeps: 300
  });

  // Both reproduce the conditioning data exactly
  for (let k = 0; k < dataIndices.length; k++) {
    assert.equal(sim[dataIndices[k]], dataValues[k]);
    assert.equal(est[dataIndices[k]], dataValues[k]);
  }

  // With noise, field fluctuates around the estimate away from data
  let diffSq = 0;
  for (let i = 0; i < w * h; i++) {
    const d = sim[i] - est[i];
    diffSq += d * d;
  }
  assert.ok(diffSq > 1.0, "Simulation with noise differs from smooth estimate away from data");

  // Verify discrete SPDE equation (kappa^2 I + L) sim == noise at free interior cells
  const lap = laplacian2D(sim, w, h, { periodic: false });
  const fixed = new Set(dataIndices);
  for (let i = 0; i < w * h; i++) {
    if (fixed.has(i)) continue;
    const lhs = 0.2 * 0.2 * sim[i] + lap[i];
    assert.ok(Math.abs(lhs - noise[i]) < 1e-4, `SPDE residual too large at cell ${i}: ${lhs} vs ${noise[i]}`);
  }
});

test("Walker Lake loader parses CSV and builds normal-score block case", () => {
  const sampleText = fs.readFileSync(new URL("../data/walker-lake/sample.csv", import.meta.url), "utf8");
  const exhText = fs.readFileSync(new URL("../data/walker-lake/exhaustive.csv", import.meta.url), "utf8");
  const sample = parseCsv(sampleText);
  const exhaustive = parseCsv(exhText);
  assert.equal(sample.V.length, 470);
  assert.equal(exhaustive.V.length, 78000);

  const c = makeCase({ sample, exhaustive }, { n: 64 });
  assert.equal(c.w, 33);
  assert.equal(c.h, 38);
  assert.ok(c.cells.length > 100);
  assert.equal(c.cells.length, c.values.length);
});
