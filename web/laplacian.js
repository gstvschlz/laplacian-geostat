// Pure numerics for the discrete graph Laplacian, neighbour relaxation,
// SPDE noise filtering, and conditional simulation in 1-D and 2-D.

/**
 * Rational approximation to the standard normal inverse CDF (Acklam's method).
 * @param {number} p - probability in (0, 1)
 * @returns {number}
 */
export function normInv(p) {
  const q = Math.min(Math.max(p, 1e-12), 1 - 1e-12);
  const a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02,
        a3 = -2.759285104469687e+02, a4 = 1.383577518672690e+02,
        a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00;
  const b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02,
        b3 = -1.556989798598866e+02, b4 = 6.680131188771972e+01,
        b5 = -1.328068155288572e+01;
  const c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01,
        c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00,
        c5 =  4.374664141464968e+00, c6 =  2.938163982698783e+00;
  const d1 =  7.784695709041462e-03, d2 =  3.224671290700398e-01,
        d3 =  2.445134137142996e+00, d4 =  3.754408661907416e+00;
  const pLow = 0.02425, pHigh = 1 - pLow;
  if (q < pLow) {
    const r = Math.sqrt(-2 * Math.log(q));
    return (((((c1 * r + c2) * r + c3) * r + c4) * r + c5) * r + c6) /
           ((((d1 * r + d2) * r + d3) * r + d4) * r + 1);
  }
  if (q <= pHigh) {
    const r = q - 0.5, s = r * r;
    return (((((a1 * s + a2) * s + a3) * s + a4) * s + a5) * s + a6) * r /
           (((((b1 * s + b2) * s + b3) * s + b4) * s + b5) * s + 1);
  }
  const r = Math.sqrt(-2 * Math.log(1 - q));
  return -(((((c1 * r + c2) * r + c3) * r + c4) * r + c5) * r + c6) /
          ((((d1 * r + d2) * r + d3) * r + d4) * r + 1);
}

/**
 * Seeded PRNG (Mulberry32) with Box-Muller normal generator.
 * @param {number} seed
 */
export function mulberry32(seed = 1) {
  let a = seed >>> 0;
  const uniform = () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  return {
    uniform,
    normal() {
      const u1 = Math.max(uniform(), 1e-12);
      const u2 = uniform();
      return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    },
    normalArray(n) {
      const out = new Float64Array(n);
      for (let i = 0; i < n; i++) out[i] = this.normal();
      return out;
    }
  };
}

/**
 * 1-D graph Laplacian: (L z)_i = d_i * z_i - sum_{j in N(i)} z_j.
 * Positive on a peak, negative in a pit, zero on a flat or a straight line.
 * @param {ArrayLike<number>} z
 * @param {{ periodic?: boolean }} [opts]
 * @returns {Float64Array}
 */
export function laplacian1D(z, { periodic = false } = {}) {
  const n = z.length;
  const out = new Float64Array(n);
  if (n < 2) return out;
  for (let i = 1; i < n - 1; i++) {
    out[i] = 2 * z[i] - z[i - 1] - z[i + 1];
  }
  if (periodic) {
    out[0] = 2 * z[0] - z[n - 1] - z[1];
    out[n - 1] = 2 * z[n - 1] - z[n - 2] - z[0];
  } else {
    out[0] = z[0] - z[1];
    out[n - 1] = z[n - 1] - z[n - 2];
  }
  return out;
}

/**
 * 1-D neighbour discrepancy: value minus average of neighbours (z_i - mean(N(i))).
 * At interior cells this is 0.5 * (L z)_i.
 * @param {ArrayLike<number>} z
 * @param {{ periodic?: boolean }} [opts]
 * @returns {Float64Array}
 */
export function discrepancy1D(z, { periodic = false } = {}) {
  const n = z.length;
  const out = new Float64Array(n);
  if (n < 2) return out;
  for (let i = 1; i < n - 1; i++) {
    out[i] = z[i] - 0.5 * (z[i - 1] + z[i + 1]);
  }
  if (periodic) {
    out[0] = z[0] - 0.5 * (z[n - 1] + z[1]);
    out[n - 1] = z[n - 1] - 0.5 * (z[n - 2] + z[0]);
  } else {
    out[0] = z[0] - z[1];
    out[n - 1] = z[n - 1] - z[n - 2];
  }
  return out;
}

/**
 * 2-D graph Laplacian: (L z)_{y,x} = d * z_{y,x} - sum of 4 neighbours.
 * @param {ArrayLike<number>} z - flattened row-major array of length width * height
 * @param {number} width
 * @param {number} height
 * @param {{ periodic?: boolean }} [opts]
 * @returns {Float64Array}
 */
export function laplacian2D(z, width, height, { periodic = false } = {}) {
  const out = new Float64Array(width * height);
  for (let y = 0; y < height; y++) {
    const ym = y > 0 ? y - 1 : (periodic ? height - 1 : -1);
    const yp = y < height - 1 ? y + 1 : (periodic ? 0 : -1);
    for (let x = 0; x < width; x++) {
      const xm = x > 0 ? x - 1 : (periodic ? width - 1 : -1);
      const xp = x < width - 1 ? x + 1 : (periodic ? 0 : -1);
      const idx = y * width + x;
      let sum = 0, deg = 0;
      if (ym !== -1) { sum += z[ym * width + x]; deg++; }
      if (yp !== -1) { sum += z[yp * width + x]; deg++; }
      if (xm !== -1) { sum += z[y * width + xm]; deg++; }
      if (xp !== -1) { sum += z[y * width + xp]; deg++; }
      out[idx] = deg * z[idx] - sum;
    }
  }
  return out;
}

/**
 * Red-black Gauss-Seidel relaxation in 1-D for (kappa^2 I + L) z = source,
 * holding dataIndices fixed at dataValues.
 * @param {Float64Array} z - state array (modified in place and returned)
 * @param {{
 *   kappa?: number,
 *   dataIndices?: ArrayLike<number>,
 *   dataValues?: ArrayLike<number>,
 *   sweeps?: number,
 *   periodic?: boolean,
 *   source?: ArrayLike<number> | null,
 *   mean?: number
 * }} [opts]
 * @returns {Float64Array}
 */
export function relax1D(z, {
  kappa = 0,
  dataIndices = [],
  dataValues = [],
  sweeps = 1,
  periodic = false,
  source = null,
  mean = 0
} = {}) {
  const n = z.length;
  const fixed = new Uint8Array(n);
  for (let k = 0; k < dataIndices.length; k++) {
    const idx = dataIndices[k];
    if (idx >= 0 && idx < n) {
      fixed[idx] = 1;
      z[idx] = dataValues[k];
    }
  }
  const k2 = kappa * kappa;
  for (let s = 0; s < sweeps; s++) {
    for (let parity = 0; parity < 2; parity++) {
      for (let i = parity; i < n; i += 2) {
        if (fixed[i]) continue;
        let nbSum = 0, deg = 0;
        if (i > 0) { nbSum += z[i - 1]; deg++; }
        else if (periodic) { nbSum += z[n - 1]; deg++; }
        if (i < n - 1) { nbSum += z[i + 1]; deg++; }
        else if (periodic) { nbSum += z[0]; deg++; }
        const src = source ? source[i] : 0;
        z[i] = (nbSum + k2 * mean + src) / (deg + k2);
      }
    }
  }
  return z;
}

/**
 * Exact solver for 1-D (kappa^2 I + L) z = source with fixed data points
 * using Thomas algorithm on the non-periodic tridiagonal system.
 * @param {number} n
 * @param {{
 *   kappa?: number,
 *   dataIndices?: ArrayLike<number>,
 *   dataValues?: ArrayLike<number>,
 *   source?: ArrayLike<number> | null,
 *   mean?: number
 * }} [opts]
 * @returns {Float64Array}
 */
export function solve1D(n, {
  kappa = 0,
  dataIndices = [],
  dataValues = [],
  source = null,
  mean = 0
} = {}) {
  const k2 = kappa * kappa;
  const a = new Float64Array(n); // sub-diagonal
  const b = new Float64Array(n); // main diagonal
  const c = new Float64Array(n); // super-diagonal
  const d = new Float64Array(n); // rhs

  for (let i = 0; i < n; i++) {
    const deg = (i > 0 ? 1 : 0) + (i < n - 1 ? 1 : 0);
    a[i] = i > 0 ? -1 : 0;
    c[i] = i < n - 1 ? -1 : 0;
    b[i] = deg + k2;
    d[i] = k2 * mean + (source ? source[i] : 0);
  }
  // If kappa == 0 and no data points are given, avoid singular matrix
  for (let k = 0; k < dataIndices.length; k++) {
    const idx = dataIndices[k];
    if (idx >= 0 && idx < n) {
      a[idx] = 0;
      c[idx] = 0;
      b[idx] = 1;
      d[idx] = dataValues[k];
    }
  }

  // Forward elimination
  for (let i = 1; i < n; i++) {
    if (b[i - 1] === 0) continue;
    const w = a[i] / b[i - 1];
    b[i] -= w * c[i - 1];
    d[i] -= w * d[i - 1];
  }
  // Back substitution
  const z = new Float64Array(n);
  z[n - 1] = b[n - 1] !== 0 ? d[n - 1] / b[n - 1] : 0;
  for (let i = n - 2; i >= 0; i--) {
    z[i] = b[i] !== 0 ? (d[i] - c[i] * z[i + 1]) / b[i] : 0;
  }
  return z;
}

/**
 * Red-black Gauss-Seidel relaxation in 2-D for (kappa^2 I + L) z = source,
 * holding dataIndices fixed at dataValues.
 * @param {Float64Array} z - flattened row-major array of length width * height
 * @param {number} width
 * @param {number} height
 * @param {{
 *   kappa?: number,
 *   dataIndices?: ArrayLike<number>,
 *   dataValues?: ArrayLike<number>,
 *   sweeps?: number,
 *   periodic?: boolean,
 *   source?: ArrayLike<number> | null,
 *   noiseScale?: number,
 *   mean?: number,
 *   omega?: number
 * }} [opts]
 * @returns {Float64Array}
 */
export function relax2D(z, width, height, {
  kappa = 0,
  dataIndices = [],
  dataValues = [],
  sweeps = 1,
  periodic = false,
  source = null,
  noiseScale = 1,
  mean = 0,
  omega = 1.0
} = {}) {
  const n = width * height;
  const fixed = new Uint8Array(n);
  for (let k = 0; k < dataIndices.length; k++) {
    const idx = dataIndices[k];
    if (idx >= 0 && idx < n) {
      fixed[idx] = 1;
      z[idx] = dataValues[k];
    }
  }
  const k2 = kappa * kappa;
  for (let s = 0; s < sweeps; s++) {
    for (let colour = 0; colour < 2; colour++) {
      for (let y = 0; y < height; y++) {
        const ym = y > 0 ? y - 1 : (periodic ? height - 1 : -1);
        const yp = y < height - 1 ? y + 1 : (periodic ? 0 : -1);
        const xStart = (y + colour) & 1;
        for (let x = xStart; x < width; x += 2) {
          const idx = y * width + x;
          if (fixed[idx]) continue;
          const xm = x > 0 ? x - 1 : (periodic ? width - 1 : -1);
          const xp = x < width - 1 ? x + 1 : (periodic ? 0 : -1);
          let nb = 0, deg = 0;
          if (ym !== -1) { nb += z[ym * width + x]; deg++; }
          if (yp !== -1) { nb += z[yp * width + x]; deg++; }
          if (xm !== -1) { nb += z[y * width + xm]; deg++; }
          if (xp !== -1) { nb += z[y * width + xp]; deg++; }
          const src = source ? noiseScale * source[idx] : 0;
          const target = (nb + k2 * mean + src) / (deg + k2);
          z[idx] = omega === 1.0 ? target : z[idx] + omega * (target - z[idx]);
        }
      }
    }
  }
  return z;
}

/**
 * Conditional simulation on a 2-D grid by solving (kappa^2 I + L) z = noise
 * with fixed data values at dataIndices. When noiseScale = 0, returns the
 * conditional mean (kriging-like estimate).
 * @param {number} width
 * @param {number} height
 * @param {{
 *   kappa?: number,
 *   dataIndices?: ArrayLike<number>,
 *   dataValues?: ArrayLike<number>,
 *   source?: ArrayLike<number> | null,
 *   noiseScale?: number,
 *   sweeps?: number,
 *   periodic?: boolean
 * }} [opts]
 * @returns {Float64Array}
 */
export function conditionalSimulate2D(width, height, {
  kappa = 0.15,
  dataIndices = [],
  dataValues = [],
  source = null,
  noiseScale = 1,
  sweeps = 250,
  periodic = false
} = {}) {
  const z = new Float64Array(width * height);
  return relax2D(z, width, height, {
    kappa,
    dataIndices,
    dataValues,
    sweeps,
    periodic,
    source,
    noiseScale
  });
}
