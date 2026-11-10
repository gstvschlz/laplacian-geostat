// Interface text written by app.js, in the language of the page (<html lang>).
// Portuguese at the root, English under en/, as on gstvschlz.com.

const pt = typeof document !== "undefined" && document.documentElement.lang.startsWith("pt");

/** Number with `d` decimals; decimal comma in Portuguese. */
export const fmt = (x, d) => {
  const s = x.toFixed(d);
  return pt ? s.replace(".", ",") : s;
};
const signed = (x, d) => (x >= 0 ? "+" : "") + fmt(x, d);

const EN = {
  value: "value",
  valMinusAvg: "val − avg",
  estimate: "estimate",
  laplacian: "Laplacian",
  mean0: "mean = 0",
  flat: "flat / straight slope (≈ 0)",
  peak: "peak (+)",
  pit: "pit (−)",
  step1: (i, v, m, d, tag) => `cell ${i}: value ${fmt(v, 2)} | neighbour avg ${fmt(m, 2)} | diff ${signed(d, 2)} → ${tag}`,
  step2: (s, h, l) => `sweeps: ${s} | held data: ${h} | max free |Laplacian|: ${fmt(l, 4)}`,
  play: "Play",
  pause: "Pause",
  rangeInf: "∞ (no sill, linear)",
  rangeCells: (r) => `≈ ${fmt(r, 1)} cells (3/κ)`,
  step3: (f, r) => `rule: sum of 2 neighbours / (2 + κ²) = ${f} × avg | practical range ${r}`,
  lapMap: (s) => `Laplacian Lz (sweep ${s})`,
  relaxed: (s) => `Relaxed estimate (sweep ${s})`,
  step4: (s, n, e) => `sweeps: ${s} | ${n} occupied blocks | normal-score RMSE vs truth: ${fmt(e, 2)}`,
  showingLap: "Showing: Laplacian Lz",
  showLap: "Show Laplacian map",
  condSim: (seed) => `Conditional simulation (seed #${seed})`,
  noiseOffLabel: "Noise off → collapses to kriging estimate",
  step6: (n, on, s) => `data held: ${n} blocks | noise: ${on ? "ON" : "OFF"} | sweeps: ${s}`,
  noiseOn: "Noise: ON",
  noiseOff: "Noise: OFF (kriging)",
  step5: (seed, s, r) => `seed #${seed} | sweeps: ${s} | Matérn range ≈ 4/κ = ${r} cells`,
  loadFailed: "Could not load the Walker Lake data.",
};

const PT = {
  value: "valor",
  valMinusAvg: "valor − média",
  estimate: "estimativa",
  laplacian: "laplaciano",
  mean0: "média = 0",
  flat: "plano / rampa reta (≈ 0)",
  peak: "pico (+)",
  pit: "vale (−)",
  step1: (i, v, m, d, tag) => `célula ${i}: valor ${fmt(v, 2)} | média dos vizinhos ${fmt(m, 2)} | diferença ${signed(d, 2)} → ${tag}`,
  step2: (s, h, l) => `passadas: ${s} | dados fixos: ${h} | máx. |laplaciano| livre: ${fmt(l, 4)}`,
  play: "Rodar",
  pause: "Pausar",
  rangeInf: "∞ (sem patamar, linear)",
  rangeCells: (r) => `≈ ${fmt(r, 1)} células (3/κ)`,
  step3: (f, r) => `regra: soma dos 2 vizinhos / (2 + κ²) = ${f} × média | alcance prático ${r}`,
  lapMap: (s) => `Laplaciano Lz (passada ${s})`,
  relaxed: (s) => `Estimativa relaxada (passada ${s})`,
  step4: (s, n, e) => `passadas: ${s} | ${n} blocos com dado | RMSE em escore normal contra a referência: ${fmt(e, 2)}`,
  showingLap: "Mostrando: laplaciano Lz",
  showLap: "Mostrar o laplaciano",
  condSim: (seed) => `Simulação condicional (semente #${seed})`,
  noiseOffLabel: "Sem ruído → volta à estimativa de krigagem",
  step6: (n, on, s) => `dados fixos: ${n} blocos | ruído: ${on ? "LIGADO" : "DESLIGADO"} | passadas: ${s}`,
  noiseOn: "Ruído: ligado",
  noiseOff: "Ruído: desligado (krigagem)",
  step5: (seed, s, r) => `semente #${seed} | passadas: ${s} | alcance Matérn ≈ 4/κ = ${r} células`,
  loadFailed: "Não foi possível carregar os dados de Walker Lake.",
};

export const t = pt ? PT : EN;
