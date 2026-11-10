import {
  discrepancy1D,
  laplacian2D,
  relax1D,
  solve1D,
  relax2D,
  mulberry32
} from "./laplacian.js";
import { divergingHex, renderGridCanvas } from "./colormap.js";
import { load as loadWalker, makeCase } from "./walker.js";

function themeColors() {
  const s = getComputedStyle(document.documentElement);
  return {
    paper: s.getPropertyValue("--color-paper").trim() || "#fcfcfb",
    ink: s.getPropertyValue("--color-ink").trim() || "#0b0b0b",
    muted: s.getPropertyValue("--color-muted").trim() || "#52514e",
    grid: s.getPropertyValue("--color-grid").trim() || "#e4e3df",
    blue: s.getPropertyValue("--color-blue").trim() || "#2a78d6",
    orange: s.getPropertyValue("--color-orange").trim() || "#eb6834"
  };
}

function fitHiDPI(canvas) {
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  const w = Math.max(200, Math.round(rect.width * dpr));
  const h = Math.max(150, Math.round(rect.height * dpr));
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, width: rect.width, height: rect.height };
}

// ---------------------------------------------------------------------------
// Step 1: The operator (1-D interactive profile + Laplacian bars)
// ---------------------------------------------------------------------------
function initStep1() {
  const canvas = document.getElementById("canvas-step1");
  const readout = document.getElementById("readout-step1");
  const N = 25;
  const z = new Float64Array(N);
  let activeIdx = 8;

  function setPreset(type) {
    for (let i = 0; i < N; i++) {
      if (type === "flat") {
        z[i] = 0.4;
      } else if (type === "slope") {
        z[i] = -1.1 + (2.2 * i) / (N - 1);
      } else if (type === "peak") {
        const p = Math.exp(-(((i - 8) / 2.2) ** 2)) * 1.45;
        const pit = -Math.exp(-(((i - 17) / 2.4) ** 2)) * 1.35;
        z[i] = -0.3 + 0.025 * i + p + pit;
      }
    }
    draw();
  }

  function xToCell(px, w) {
    const padL = 42, padR = 16;
    const plotW = w - padL - padR;
    const idx = Math.round(((px - padL) / plotW) * (N - 1));
    return Math.max(0, Math.min(N - 1, idx));
  }

  function yToVal(py, h) {
    const topY = 18, topH = h * 0.54;
    const u = 1 - (py - topY) / topH;
    const val = -1.8 + u * 3.6;
    return Math.max(-1.7, Math.min(1.7, val));
  }

  function draw() {
    const { ctx, width: w, height: h } = fitHiDPI(canvas);
    const c = themeColors();
    ctx.clearRect(0, 0, w, h);

    const padL = 42, padR = 16;
    const plotW = w - padL - padR;
    const topY = 18, topH = h * 0.52;
    const botY = topY + topH + 22, botH = h - botY - 18;

    const cellX = (i) => padL + (i / (N - 1)) * plotW;
    const valY = (v) => topY + (1 - (v + 1.8) / 3.6) * topH;
    const lapZeroY = botY + botH * 0.5;

    // Zero lines & labels
    ctx.strokeStyle = c.grid;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, valY(0));
    ctx.lineTo(w - padR, valY(0));
    ctx.moveTo(padL, lapZeroY);
    ctx.lineTo(w - padR, lapZeroY);
    ctx.stroke();

    ctx.fillStyle = c.muted;
    ctx.font = "italic 14px \"Crimson Pro\", Georgia, serif";
    ctx.textAlign = "right";
    ctx.textAlign = "left"; ctx.fillText("value", 2, topY - 4); ctx.textAlign = "right";
    ctx.fillText("0", padL - 6, valY(0) + 4);
    ctx.textAlign = "left"; ctx.fillText("val − avg", 2, botY - 4); ctx.textAlign = "right";
    ctx.fillText("0", padL - 6, lapZeroY + 4);

    const disc = discrepancy1D(z);

    // Profile polyline
    ctx.strokeStyle = c.grid;
    ctx.lineWidth = 1.75;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = cellX(i), y = valY(z[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Highlight active interior cell's neighbour chord and difference bar
    const hi = Math.max(1, Math.min(N - 2, activeIdx));
    const nbMean = 0.5 * (z[hi - 1] + z[hi + 1]);
    ctx.strokeStyle = c.blue;
    ctx.lineWidth = 2;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(cellX(hi - 1), valY(z[hi - 1]));
    ctx.lineTo(cellX(hi + 1), valY(z[hi + 1]));
    ctx.stroke();
    ctx.setLineDash([]);

    // Vertical difference bar on profile
    ctx.strokeStyle = disc[hi] >= 0 ? c.orange : c.blue;
    ctx.lineWidth = 3.5;
    ctx.beginPath();
    ctx.moveTo(cellX(hi), valY(nbMean));
    ctx.lineTo(cellX(hi), valY(z[hi]));
    ctx.stroke();

    // Midpoint hollow circle
    ctx.fillStyle = c.paper;
    ctx.strokeStyle = c.blue;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(cellX(hi), valY(nbMean), 4.5, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();

    // Draw nodes and bottom Laplacian bars
    const barW = Math.max(4, (plotW / N) * 0.68);
    for (let i = 0; i < N; i++) {
      const x = cellX(i);
      const y = valY(z[i]);
      const d = disc[i];
      const col = Math.abs(d) < 0.015 ? c.muted : divergingHex(d / 0.9);

      // Bottom bar
      const bh = Math.max(-botH * 0.48, Math.min(botH * 0.48, (d / 1.2) * (botH * 0.48)));
      ctx.fillStyle = Math.abs(d) < 0.015 ? c.grid : (d > 0 ? c.orange : c.blue);
      if (bh >= 0) {
        ctx.fillRect(x - barW / 2, lapZeroY - bh, barW, Math.max(1.5, bh));
      } else {
        ctx.fillRect(x - barW / 2, lapZeroY, barW, Math.max(1.5, -bh));
      }

      // Top node
      ctx.fillStyle = col;
      ctx.strokeStyle = i === hi ? c.ink : c.paper;
      ctx.lineWidth = i === hi ? 2 : 1.2;
      ctx.beginPath();
      ctx.arc(x, y, i === hi ? 6.5 : 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }

    const diff = z[hi] - nbMean;
    const tag = Math.abs(diff) < 0.02 ? "flat / straight slope (≈ 0)" : (diff > 0 ? "peak (+)" : "pit (−)");
    readout.textContent = `cell ${hi}: value ${z[hi].toFixed(2)} | neighbour avg ${nbMean.toFixed(2)} | diff ${diff >= 0 ? "+" : ""}${diff.toFixed(2)} → ${tag}`;
  }

  let dragging = false;
  canvas.addEventListener("pointerdown", (e) => {
    dragging = true;
    canvas.setPointerCapture(e.pointerId);
    const rect = canvas.getBoundingClientRect();
    const idx = xToCell(e.clientX - rect.left, rect.width);
    activeIdx = idx;
    z[idx] = yToVal(e.clientY - rect.top, rect.height);
    draw();
  });
  canvas.addEventListener("pointermove", (e) => {
    const rect = canvas.getBoundingClientRect();
    const idx = xToCell(e.clientX - rect.left, rect.width);
    if (!dragging) {
      if (idx !== activeIdx) {
        activeIdx = idx;
        draw();
      }
      return;
    }
    activeIdx = idx;
    z[idx] = yToVal(e.clientY - rect.top, rect.height);
    draw();
  });
  canvas.addEventListener("pointerup", () => { dragging = false; });

  document.querySelectorAll("[data-step1-preset]").forEach((btn) => {
    btn.addEventListener("click", () => setPreset(btn.dataset.step1Preset));
  });

  setPreset("peak");
  return { draw };
}

// ---------------------------------------------------------------------------
// Step 2: Estimation by repeated neighbour averaging (1-D)
// ---------------------------------------------------------------------------
function initStep2() {
  const canvas = document.getElementById("canvas-step2");
  const readout = document.getElementById("readout-step2");
  const playBtn = document.getElementById("btn-step2-play");
  const stepBtn = document.getElementById("btn-step2-step");
  const resetBtn = document.getElementById("btn-step2-reset");

  const N = 41;
  const z = new Float64Array(N);
  const held = new Map([
    [4, 0.5],
    [13, 1.35],
    [21, 0.45],
    [30, -1.05],
    [37, 0.4]
  ]);
  let sweeps = 0;
  let playing = false;

  function syncHeld() {
    for (const [idx, val] of held.entries()) z[idx] = val;
  }

  function resetField() {
    z.fill(0);
    syncHeld();
    sweeps = 0;
    draw();
  }

  function stepOnce(count = 1) {
    const dataIndices = Array.from(held.keys());
    const dataValues = Array.from(held.values());
    relax1D(z, { kappa: 0, dataIndices, dataValues, sweeps: count, periodic: false });
    sweeps += count;
    draw();
  }

  function draw() {
    const { ctx, width: w, height: h } = fitHiDPI(canvas);
    const c = themeColors();
    ctx.clearRect(0, 0, w, h);

    const padL = 42, padR = 16;
    const plotW = w - padL - padR;
    const topY = 18, topH = h * 0.54;
    const botY = topY + topH + 22, botH = h - botY - 16;

    const cellX = (i) => padL + (i / (N - 1)) * plotW;
    const valY = (v) => topY + (1 - (v + 1.6) / 3.2) * topH;
    const lapZeroY = botY + botH * 0.5;

    ctx.strokeStyle = c.grid;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, valY(0));
    ctx.lineTo(w - padR, valY(0));
    ctx.moveTo(padL, lapZeroY);
    ctx.lineTo(w - padR, lapZeroY);
    ctx.stroke();

    ctx.fillStyle = c.muted;
    ctx.font = "italic 14px \"Crimson Pro\", Georgia, serif";
    ctx.textAlign = "right";
    ctx.textAlign = "left"; ctx.fillText("estimate", 2, topY - 4); ctx.textAlign = "right";
    ctx.textAlign = "left"; ctx.fillText("Laplacian", 2, botY - 4); ctx.textAlign = "right";

    // Current estimate curve
    ctx.strokeStyle = c.blue;
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = cellX(i), y = valY(z[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    const disc = discrepancy1D(z);
    const barW = Math.max(3, (plotW / N) * 0.65);
    let maxFreeLap = 0;

    for (let i = 0; i < N; i++) {
      const x = cellX(i);
      const isHeld = held.has(i);
      const d = isHeld ? 0 : disc[i];
      if (!isHeld && Math.abs(d) > maxFreeLap) maxFreeLap = Math.abs(d);

      const bh = Math.max(-botH * 0.46, Math.min(botH * 0.46, (d / 0.6) * (botH * 0.46)));
      ctx.fillStyle = d >= 0 ? c.orange : c.blue;
      if (Math.abs(d) > 0.002) {
        if (bh >= 0) ctx.fillRect(x - barW / 2, lapZeroY - bh, barW, bh);
        else ctx.fillRect(x - barW / 2, lapZeroY, barW, -bh);
      }

      ctx.beginPath();
      ctx.arc(x, valY(z[i]), isHeld ? 6.5 : 3.2, 0, Math.PI * 2);
      ctx.fillStyle = isHeld ? c.ink : c.blue;
      ctx.strokeStyle = c.paper;
      ctx.lineWidth = 1.5;
      ctx.fill();
      if (isHeld) ctx.stroke();
    }

    readout.textContent = `sweeps: ${sweeps} | held data: ${held.size} | max free |Laplacian|: ${maxFreeLap.toFixed(4)}`;
  }

  let draggingCell = null;
  canvas.addEventListener("pointerdown", (e) => {
    const rect = canvas.getBoundingClientRect();
    const padL = 42, padR = 16;
    const plotW = rect.width - padL - padR;
    const topY = 18, topH = rect.height * 0.54;
    const idx = Math.max(0, Math.min(N - 1, Math.round(((e.clientX - rect.left - padL) / plotW) * (N - 1))));
    const val = Math.max(-1.5, Math.min(1.5, -1.6 + (1 - (e.clientY - rect.top - topY) / topH) * 3.2));

    if (e.shiftKey && held.has(idx) && held.size > 1) {
      held.delete(idx);
      draw();
      return;
    }
    held.set(idx, val);
    z[idx] = val;
    draggingCell = idx;
    canvas.setPointerCapture(e.pointerId);
    draw();
  });

  canvas.addEventListener("pointermove", (e) => {
    if (draggingCell === null) return;
    const rect = canvas.getBoundingClientRect();
    const topY = 18, topH = rect.height * 0.54;
    const val = Math.max(-1.5, Math.min(1.5, -1.6 + (1 - (e.clientY - rect.top - topY) / topH) * 3.2));
    held.set(draggingCell, val);
    z[draggingCell] = val;
    draw();
  });

  canvas.addEventListener("pointerup", () => { draggingCell = null; });
  canvas.addEventListener("dblclick", (e) => {
    const rect = canvas.getBoundingClientRect();
    const padL = 42, padR = 16;
    const plotW = rect.width - padL - padR;
    const idx = Math.max(0, Math.min(N - 1, Math.round(((e.clientX - rect.left - padL) / plotW) * (N - 1))));
    if (held.has(idx) && held.size > 1) {
      held.delete(idx);
      draw();
    }
  });

  playBtn.addEventListener("click", () => {
    playing = !playing;
    playBtn.setAttribute("aria-pressed", String(playing));
    playBtn.textContent = playing ? "Pause" : "Play";
  });

  stepBtn.addEventListener("click", () => stepOnce(1));
  resetBtn.addEventListener("click", () => resetField());

  function tick() {
    if (playing) stepOnce(2);
    requestAnimationFrame(tick);
  }
  resetField();
  requestAnimationFrame(tick);
  return { draw };
}

// ---------------------------------------------------------------------------
// Step 3: kappa slider — sagging toward the mean
// ---------------------------------------------------------------------------
function initStep3() {
  const canvas = document.getElementById("canvas-step3");
  const slider = document.getElementById("slider-step3-kappa");
  const valSpan = document.getElementById("val-step3-kappa");
  const readout = document.getElementById("readout-step3");

  const N = 49;
  const held = new Map([
    [4, 0.55],
    [14, 1.35],
    [22, 0.65],
    [33, -0.95],
    [44, 0.45]
  ]);

  function draw() {
    const kappa = parseFloat(slider.value);
    valSpan.textContent = kappa.toFixed(2);

    const dataIndices = Array.from(held.keys());
    const dataValues = Array.from(held.values());
    const linear = solve1D(N, { kappa: 0, dataIndices, dataValues });
    const sagged = solve1D(N, { kappa, dataIndices, dataValues });

    const { ctx, width: w, height: h } = fitHiDPI(canvas);
    const c = themeColors();
    ctx.clearRect(0, 0, w, h);

    const padL = 42, padR = 16, topY = 20, plotH = h - 42;
    const plotW = w - padL - padR;
    const cellX = (i) => padL + (i / (N - 1)) * plotW;
    const valY = (v) => topY + (1 - (v + 1.5) / 3.1) * plotH;

    // Mean line (0)
    ctx.strokeStyle = c.muted;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(padL, valY(0));
    ctx.lineTo(w - padR, valY(0));
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = c.muted;
    ctx.font = "italic 14px \"Crimson Pro\", Georgia, serif";
    ctx.textAlign = "right";
    ctx.fillText("mean = 0", w - padR - 4, valY(0) - 6);

    // Reference kappa=0 straight-line interpolation
    ctx.strokeStyle = c.grid;
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = cellX(i), y = valY(linear[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Sagged curve for current kappa
    ctx.strokeStyle = c.blue;
    ctx.lineWidth = 2.6;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = cellX(i), y = valY(sagged[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // Held data points
    for (const [idx, val] of held.entries()) {
      ctx.beginPath();
      ctx.arc(cellX(idx), valY(val), 6.5, 0, Math.PI * 2);
      ctx.fillStyle = c.ink;
      ctx.strokeStyle = c.paper;
      ctx.lineWidth = 1.5;
      ctx.fill();
      ctx.stroke();
    }

    const factor = (2 / (2 + kappa * kappa)).toFixed(3);
    const rangeText = kappa < 0.02 ? "∞ (no sill, linear)" : `≈ ${(3 / kappa).toFixed(1)} cells (3/κ)`;
    readout.textContent = `rule: sum of 2 neighbours / (2 + κ²) = ${factor} × avg | practical range ${rangeText}`;
  }

  let draggingIdx = null;
  canvas.addEventListener("pointerdown", (e) => {
    const rect = canvas.getBoundingClientRect();
    const padL = 42, padR = 16, topY = 20, plotH = rect.height - 42;
    const plotW = rect.width - padL - padR;
    const idx = Math.max(0, Math.min(N - 1, Math.round(((e.clientX - rect.left - padL) / plotW) * (N - 1))));
    const val = Math.max(-1.4, Math.min(1.5, -1.5 + (1 - (e.clientY - rect.top - topY) / plotH) * 3.1));
    held.set(idx, val);
    draggingIdx = idx;
    canvas.setPointerCapture(e.pointerId);
    draw();
  });
  canvas.addEventListener("pointermove", (e) => {
    if (draggingIdx === null) return;
    const rect = canvas.getBoundingClientRect();
    const topY = 20, plotH = rect.height - 42;
    const val = Math.max(-1.4, Math.min(1.5, -1.5 + (1 - (e.clientY - rect.top - topY) / plotH) * 3.1));
    held.set(draggingIdx, val);
    draw();
  });
  canvas.addEventListener("pointerup", () => { draggingIdx = null; });

  slider.addEventListener("input", draw);
  draw();
  return { draw };
}

// ---------------------------------------------------------------------------
// Step 4 & Step 6: Walker Lake 2-D Map & Conditional Simulation
// ---------------------------------------------------------------------------
async function initWalkerSteps() {
  let walkerCase;
  try {
    const raw = await loadWalker("data/walker-lake");
    walkerCase = makeCase(raw, { n: 64 });
  } catch (err) {
    console.error("Failed to load Walker Lake CSVs:", err);
    return;
  }

  const N = 64;
  const { w, h, truth, cells, values } = walkerCase;

  // --- Step 4: Map estimation ---
  const canvasEst = document.getElementById("canvas-step4-est");
  const canvasTruth = document.getElementById("canvas-step4-truth");
  const sliderK4 = document.getElementById("slider-step4-kappa");
  const valK4 = document.getElementById("val-step4-kappa");
  const btnPlay4 = document.getElementById("btn-step4-play");
  const btnReset4 = document.getElementById("btn-step4-reset");
  const btnLap4 = document.getElementById("btn-step4-lap");
  const readout4 = document.getElementById("readout-step4");
  const labelLeft4 = document.getElementById("label-step4-left");

  const z4 = new Float64Array(N * N);
  let sweeps4 = 0;
  let playing4 = true;
  let showLap4 = false;

  renderGridCanvas(canvasTruth, truth, N, w, h, { limit: 2.5, samples: { cells, values } });

  function reset4() {
    z4.fill(0);
    for (let k = 0; k < cells.length; k++) z4[cells[k]] = values[k];
    sweeps4 = 0;
    draw4();
  }

  function draw4() {
    if (showLap4) {
      const lap = laplacian2D(z4, N, N, { periodic: true });
      renderGridCanvas(canvasEst, lap, N, w, h, { limit: 1.5, samples: { cells, values } });
      labelLeft4.textContent = `Laplacian Lz (sweep ${sweeps4})`;
    } else {
      renderGridCanvas(canvasEst, z4, N, w, h, { limit: 2.5, samples: { cells, values } });
      labelLeft4.textContent = `Relaxed estimate (sweep ${sweeps4})`;
    }
    // Compute RMSE against exhaustive truth in normal-score units over the w x h window
    let mse = 0, count = 0;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        const d = z4[y * N + x] - truth[y * N + x];
        mse += d * d;
        count++;
      }
    }
    const rmse = Math.sqrt(mse / count);
    readout4.textContent = `sweeps: ${sweeps4} | ${cells.length} occupied blocks | normal-score RMSE vs truth: ${rmse.toFixed(2)}`;
  }

  btnPlay4.addEventListener("click", () => {
    playing4 = !playing4;
    btnPlay4.setAttribute("aria-pressed", String(playing4));
    btnPlay4.textContent = playing4 ? "Pause" : "Play";
  });
  btnReset4.addEventListener("click", reset4);
  btnLap4.addEventListener("click", () => {
    showLap4 = !showLap4;
    btnLap4.setAttribute("aria-pressed", String(showLap4));
    btnLap4.textContent = showLap4 ? "Showing: Laplacian Lz" : "Show Laplacian map";
    draw4();
  });
  sliderK4.addEventListener("input", () => {
    valK4.textContent = parseFloat(sliderK4.value).toFixed(2);
    playing4 = true;
    btnPlay4.setAttribute("aria-pressed", "true");
    btnPlay4.textContent = "Pause";
  });

  // --- Step 6: Conditional simulation ---
  const canvasCond = document.getElementById("canvas-step6-sim");
  const canvasKrig = document.getElementById("canvas-step6-krig");
  const btnNoise6 = document.getElementById("btn-step6-noise");
  const btnSeed6 = document.getElementById("btn-step6-seed");
  const sliderK6 = document.getElementById("slider-step6-kappa");
  const valK6 = document.getElementById("val-step6-kappa");
  const readout6 = document.getElementById("readout-step6");
  const labelRight6 = document.getElementById("label-step6-right");

  let seed6 = 101;
  let noiseOn6 = true;
  let noise6 = mulberry32(seed6).normalArray(N * N);
  const zCond = new Float64Array(N * N);
  const zKrig = new Float64Array(N * N);
  let sweeps6 = 0;

  function recomputeStep6(initialSweeps = 60) {
    const kappa = parseFloat(sliderK6.value);
    valK6.textContent = kappa.toFixed(2);
    // Scale white noise so the stationary variance of (kappa^2 I + L) z = w is ~1
    const amp = Math.sqrt(4 * Math.PI) * Math.max(kappa, 0.08) * 1.35;
    relax2D(zKrig, N, N, {
      kappa,
      dataIndices: cells,
      dataValues: values,
      sweeps: initialSweeps,
      periodic: true,
      source: null
    });
    relax2D(zCond, N, N, {
      kappa,
      dataIndices: cells,
      dataValues: values,
      sweeps: initialSweeps,
      periodic: true,
      source: noise6,
      noiseScale: noiseOn6 ? amp : 0
    });
    sweeps6 += initialSweeps;
    draw6();
  }

  function draw6() {
    renderGridCanvas(canvasKrig, zKrig, N, w, h, { limit: 2.5, samples: { cells, values } });
    renderGridCanvas(canvasCond, zCond, N, w, h, { limit: 2.5, samples: { cells, values } });
    labelRight6.textContent = noiseOn6
      ? `Conditional simulation (seed #${seed6})`
      : "Noise off → collapses to kriging estimate";
    readout6.textContent = `data held: ${cells.length} blocks | noise: ${noiseOn6 ? "ON" : "OFF"} | sweeps: ${sweeps6}`;
  }

  btnNoise6.addEventListener("click", () => {
    noiseOn6 = !noiseOn6;
    btnNoise6.setAttribute("aria-pressed", String(noiseOn6));
    btnNoise6.textContent = noiseOn6 ? "Noise: ON" : "Noise: OFF (Kriging)";
    sweeps6 = 0;
    recomputeStep6(80);
  });

  btnSeed6.addEventListener("click", () => {
    seed6 += 1;
    noise6 = mulberry32(seed6).normalArray(N * N);
    noiseOn6 = true;
    btnNoise6.setAttribute("aria-pressed", "true");
    btnNoise6.textContent = "Noise: ON";
    sweeps6 = 0;
    recomputeStep6(80);
  });

  sliderK6.addEventListener("input", () => {
    sweeps6 = 0;
    recomputeStep6(60);
  });

  reset4();
  recomputeStep6(120);

  function tick() {
    if (playing4 && sweeps4 < 600) {
      const kappa4 = parseFloat(sliderK4.value);
      relax2D(z4, N, N, {
        kappa: kappa4,
        dataIndices: cells,
        dataValues: values,
        sweeps: 3,
        periodic: true
      });
      sweeps4 += 3;
      draw4();
    }
    if (sweeps6 < 360) {
      recomputeStep6(4);
    }
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}

// ---------------------------------------------------------------------------
// Step 5: Unconditional Simulation from white noise: (kappa^2 I + L) z = noise
// ---------------------------------------------------------------------------
function initStep5() {
  const canvasNoise = document.getElementById("canvas-step5-noise");
  const canvasSim = document.getElementById("canvas-step5-sim");
  const btnSeed = document.getElementById("btn-step5-seed");
  const btnPlay = document.getElementById("btn-step5-play");
  const sliderK = document.getElementById("slider-step5-kappa");
  const valK = document.getElementById("val-step5-kappa");
  const readout = document.getElementById("readout-step5");

  const N = 48;
  let seed = 7;
  let noise = mulberry32(seed).normalArray(N * N);
  const z = new Float64Array(N * N);
  let sweeps = 0;
  let playing = true;

  function renderNormalized(canvas, arr) {
    let mean = 0;
    for (let i = 0; i < arr.length; i++) mean += arr[i];
    mean /= arr.length;
    let v = 0;
    for (let i = 0; i < arr.length; i++) {
      const d = arr[i] - mean;
      v += d * d;
    }
    const std = Math.sqrt(v / arr.length) || 1;
    const norm = Float64Array.from(arr, (x) => (x - mean) / std);
    renderGridCanvas(canvas, norm, N, N, N, { limit: 2.4 });
  }

  function restart() {
    z.fill(0);
    sweeps = 0;
    renderNormalized(canvasNoise, noise);
    stepSim(2);
  }

  function stepSim(count) {
    const kappa = parseFloat(sliderK.value);
    valK.textContent = kappa.toFixed(2);
    relax2D(z, N, N, {
      kappa,
      sweeps: count,
      periodic: true,
      source: noise,
      noiseScale: 1
    });
    sweeps += count;
    renderNormalized(canvasSim, z);
    const rangeCells = Math.round(4 / kappa);
    readout.textContent = `seed #${seed} | sweeps: ${sweeps} | Matérn range ≈ 4/κ = ${rangeCells} cells`;
  }

  btnSeed.addEventListener("click", () => {
    seed += 1;
    noise = mulberry32(seed).normalArray(N * N);
    playing = true;
    btnPlay.setAttribute("aria-pressed", "true");
    btnPlay.textContent = "Pause";
    restart();
  });

  btnPlay.addEventListener("click", () => {
    playing = !playing;
    btnPlay.setAttribute("aria-pressed", String(playing));
    btnPlay.textContent = playing ? "Pause" : "Play";
  });

  sliderK.addEventListener("input", () => {
    playing = true;
    btnPlay.setAttribute("aria-pressed", "true");
    btnPlay.textContent = "Pause";
    stepSim(6);
  });

  restart();
  function tick() {
    if (playing && sweeps < 450) {
      stepSim(3);
    }
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}

// ---------------------------------------------------------------------------
// Boot all 6 steps + theme toggle
// ---------------------------------------------------------------------------
const s1 = initStep1();
const s2 = initStep2();
const s3 = initStep3();
initStep5();
initWalkerSteps();

// Theme buttons, same markup and storage key as gstvschlz.com.
const root = document.documentElement;
const themeButtons = [...document.querySelectorAll("[data-theme-set]")];
const syncThemeButtons = () =>
  themeButtons.forEach((b) => b.setAttribute("aria-pressed", String((root.dataset.theme || "") === b.dataset.themeSet)));
const redraw = () => { s1.draw(); s2.draw(); s3.draw(); };
themeButtons.forEach((b) =>
  b.addEventListener("click", () => {
    const v = b.dataset.themeSet;
    if (v) root.dataset.theme = v; else delete root.dataset.theme;
    try { v ? localStorage.setItem("theme", v) : localStorage.removeItem("theme"); } catch {}
    syncThemeButtons();
    redraw();
  })
);
syncThemeButtons();
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);

window.addEventListener("resize", () => {
  s1.draw();
  s2.draw();
  s3.draw();
});
