"""Render the README and the chapters into one page, _site/index.html."""

import re
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parents[1]
CSS = """
:root { color-scheme: light dark; --bg: #fcfcfb; --ink: #0b0b0b; --muted: #52514e; --line: #e4e3df; --code: #f1f0ec; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #161615; --ink: #f0efec; --muted: #c3c2b7; --line: #33322f; --code: #232321; }
}
body { background: var(--bg); color: var(--ink); font: 17px/1.6 Georgia, serif; margin: 0; padding: 0 16px; }
main { max-width: 46rem; margin: 0 auto; padding: 3rem 0 6rem; }
h1, h2, h3 { font-family: system-ui, sans-serif; line-height: 1.25; }
h1 { font-size: 1.9rem; margin-top: 4rem; padding-top: 2rem; border-top: 1px solid var(--line); }
section:first-child h1 { margin-top: 0; padding-top: 0; border-top: 0; }
h2 { font-size: 1.2rem; margin-top: 2.2rem; }
a { color: inherit; }
img { display: block; width: min(100vw - 32px, 62rem); max-width: none; margin: 1.5rem 0 1.5rem 50%;
      transform: translateX(-50%); border-radius: 6px; }
code { font: 0.86em ui-monospace, monospace; background: var(--code); padding: 0.1em 0.3em; border-radius: 4px; }
pre { background: var(--code); padding: 0.9rem 1rem; border-radius: 6px; overflow-x: auto; line-height: 1.45; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; font-size: 0.93em; display: block; overflow-x: auto; }
th, td { border-bottom: 1px solid var(--line); padding: 0.4rem 0.8rem 0.4rem 0; text-align: left; vertical-align: top; }
th { font-family: system-ui, sans-serif; color: var(--muted); font-weight: 600; }
"""


def render(path, anchor):
    text = path.read_text()
    folder = path.parent.relative_to(ROOT).as_posix()
    text = re.sub(r"\]\((?!https?:|#)([^)]+\.png)\)", rf"](../{folder}/\1)", text)
    text = re.sub(r"\]\((?:\.\./|chapters/)(\d\d-[a-z-]+)/\)", r"](#\1)", text)
    text = text.replace("](../../README.md)", "](#overview)")
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    return f'<section id="{anchor}">{body}</section>'


def main():
    parts = [render(ROOT / "README.md", "overview")]
    parts += [render(p, p.parent.name) for p in sorted(ROOT.glob("chapters/*/README.md"))]
    out = ROOT / "_site" / "index.html"
    out.parent.mkdir(exist_ok=True)
    out.write_text(
        '<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>The Laplacian in geostatistics</title><style>{CSS}</style><main>{''.join(parts)}</main>"
    )
    print(out.relative_to(ROOT))


if __name__ == "__main__":
    main()
