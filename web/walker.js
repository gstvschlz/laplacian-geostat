// The Walker Lake data set (Isaaks and Srivastava, 1989): 470 samples of V
// and the exhaustive 260 x 300 m grid they were drawn from. The grid is
// averaged into square blocks and placed in a corner of the torus; the rest
// of the torus is free padding, so the wrap-around does not tie opposite
// edges together. V is skewed, so the automaton works on normal scores from
// the declustered sample distribution.

import { normInv } from "./laplacian.js";

export const WIDTH = 260, HEIGHT = 300;
export const SAMPLES = 470;
export const DECLUSTER_CELL = 20;
const BLOCK = { 64: 8, 128: 4, 256: 2, 512: 1 };

export function blockSize(n) {
  const b = BLOCK[n];
  if (!b) throw new Error(`no block size for a ${n} grid`);
  return b;
}

export function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const head = lines[0].split(",");
  const cols = Object.fromEntries(head.map((h) => [h, new Float64Array(lines.length - 1)]));
  for (let r = 1; r < lines.length; r++) {
    const f = lines[r].split(",");
    for (let c = 0; c < head.length; c++) cols[head[c]][r - 1] = f[c] === "" ? NaN : Number(f[c]);
  }
  return cols;
}

export async function load(base = "data/walker-lake") {
  const get = async (name) => {
    const res = await fetch(`${base}/${name}`);
    if (!res.ok) throw new Error(`${name}: ${res.status}`);
    return parseCsv(await res.text());
  };
  const [sample, exhaustive] = await Promise.all([get("sample.csv"), get("exhaustive.csv")]);
  return { sample, exhaustive };
}

// Cell declustering: each occupied cell gets the same total weight.
export function declusterWeights(x, y, cell = DECLUSTER_CELL) {
  const keys = Array.from(x, (_, k) => `${Math.floor(x[k] / cell)},${Math.floor(y[k] / cell)}`);
  const counts = new Map();
  for (const key of keys) counts.set(key, (counts.get(key) ?? 0) + 1);
  return Float64Array.from(keys, (key) => 1 / (counts.get(key) * counts.size));
}

// Normal-score transform from a weighted sample. Tied values share the
// middle of their probability step; between sample values the cumulative
// probability is interpolated linearly, beyond them it is held.
export function anamorphosis(values, weights) {
  const order = Array.from(values, (_, k) => k).sort((a, b) => values[a] - values[b]);
  const total = order.reduce((s, k) => s + weights[k], 0);
  const v = [], p = [];
  let below = 0;
  for (let q = 0; q < order.length;) {
    const value = values[order[q]];
    let step = 0;
    for (; q < order.length && values[order[q]] === value; q++) step += weights[order[q]];
    v.push(value);
    p.push((below + step / 2) / total);
    below += step;
  }
  return (x) => {
    if (x <= v[0]) return normInv(p[0]);
    if (x >= v[v.length - 1]) return normInv(p[p.length - 1]);
    let lo = 0, hi = v.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (v[mid] <= x) lo = mid; else hi = mid; }
    return normInv(p[lo] + ((x - v[lo]) / (v[hi] - v[lo])) * (p[hi] - p[lo]));
  };
}

// Cells of an n x n torus: the data window takes columns 0..w-1 and rows
// 0..h-1, north up. Grid nodes and samples sit at integer metres from 1.
export function makeCase({ sample, exhaustive }, { n, count = SAMPLES }) {
  const b = blockSize(n);
  const w = Math.ceil(WIDTH / b), h = Math.ceil(HEIGHT / b);
  const cellOf = (x, y) => (h - 1 - Math.floor((y - 1) / b)) * n + Math.floor((x - 1) / b);
  const forward = anamorphosis(sample.V, declusterWeights(sample.X, sample.Y));

  const sum = new Float64Array(n * n), hits = new Float64Array(n * n);
  for (let k = 0; k < exhaustive.V.length; k++) {
    const c = cellOf(exhaustive.X[k], exhaustive.Y[k]);
    sum[c] += exhaustive.V[k];
    hits[c]++;
  }
  const truth = new Float64Array(n * n);
  for (let c = 0; c < n * n; c++) if (hits[c]) truth[c] = forward(sum[c] / hits[c]);

  // Samples in one block are averaged in normal-score units.
  const scores = new Map();
  for (let k = 0; k < Math.min(count, sample.V.length); k++) {
    const c = cellOf(sample.X[k], sample.Y[k]);
    const s = scores.get(c) ?? [0, 0];
    s[0] += forward(sample.V[k]);
    s[1]++;
    scores.set(c, s);
  }
  const cells = Int32Array.from(scores.keys());
  const values = Float64Array.from(scores.values(), ([s, m]) => s / m);
  return { block: b, w, h, truth, cells, values, forward };
}
