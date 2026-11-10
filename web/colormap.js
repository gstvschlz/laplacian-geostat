// Grey scale, black (low) to white (high), to match the site's black-white-grey theme.
// The names say "diverging" because u runs over [-1, 1], centred on mid grey.

const STOPS = [
  [0x10, 0x10, 0x10],
  [0x58, 0x58, 0x58],
  [0x9a, 0x9a, 0x9a],
  [0xd2, 0xd2, 0xd2],
  [0xfa, 0xfa, 0xfa]
];

/**
 * Interpolate diverging colormap for normalised value u in [-1, 1].
 * @param {number} u - value in [-1, 1]
 * @returns {[number, number, number]} [r, g, b] in 0..255
 */
export function divergingRgb(u) {
  const clamped = Math.max(-1, Math.min(1, Number.isFinite(u) ? u : 0));
  const pos = (clamped + 1) * 0.5 * (STOPS.length - 1);
  const idx = Math.min(Math.floor(pos), STOPS.length - 2);
  const frac = pos - idx;
  const a = STOPS[idx], b = STOPS[idx + 1];
  return [
    Math.round(a[0] + (b[0] - a[0]) * frac),
    Math.round(a[1] + (b[1] - a[1]) * frac),
    Math.round(a[2] + (b[2] - a[2]) * frac)
  ];
}

export function divergingHex(u) {
  const [r, g, b] = divergingRgb(u);
  return `rgb(${r}, ${g}, ${b})`;
}

/**
 * Render a 2-D scalar subgrid (0..viewW-1, 0..viewH-1) from a stride-`stride`
 * array into a canvas, plus optional sample markers.
 */
export function renderGridCanvas(canvas, field, stride, viewW, viewH, {
  limit = 2.5,
  samples = null
} = {}) {
  const scale = 6;
  const cw = viewW * scale;
  const ch = viewH * scale;
  if (canvas.width !== cw || canvas.height !== ch) {
    canvas.width = cw;
    canvas.height = ch;
  }
  const ctx = canvas.getContext("2d");
  const img = ctx.createImageData(cw, ch);
  const data = img.data;

  for (let y = 0; y < viewH; y++) {
    for (let x = 0; x < viewW; x++) {
      const val = field[y * stride + x];
      const [r, g, b] = divergingRgb(val / limit);
      for (let dy = 0; dy < scale; dy++) {
        const rowOff = ((y * scale + dy) * cw + x * scale) * 4;
        for (let dx = 0; dx < scale; dx++) {
          const p = rowOff + dx * 4;
          data[p] = r;
          data[p + 1] = g;
          data[p + 2] = b;
          data[p + 3] = 255;
        }
      }
    }
  }
  ctx.putImageData(img, 0, 0);

  if (samples && samples.cells && samples.values) {
    const { cells, values } = samples;
    ctx.lineWidth = 1.25;
    for (let k = 0; k < cells.length; k++) {
      const c = cells[k];
      const y = Math.floor(c / stride);
      const x = c % stride;
      if (x >= viewW || y >= viewH) continue;
      const cx = (x + 0.5) * scale;
      const cy = (y + 0.5) * scale;
      // A small black dot with a white rim marks a sample on any grey; the cell keeps its value.
      ctx.beginPath();
      ctx.arc(cx, cy, Math.max(1.6, scale * 0.2), 0, Math.PI * 2);
      ctx.fillStyle = "#000000";
      ctx.fill();
      ctx.lineWidth = 1;
      ctx.strokeStyle = "#ffffff";
      ctx.stroke();
    }
  }
}
