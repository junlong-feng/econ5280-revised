#!/usr/bin/env python3
"""Lightweight structural checks for the revised course sources."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
errors: list[str] = []
warnings: list[str] = []


def complain(path: Path, message: str) -> None:
    errors.append(f"{path.relative_to(ROOT)}: {message}")


tex_files = sorted((ROOT / "latex").rglob("*.tex"))
qmd_files = [ROOT / "live" / "index.qmd"]
qmd_files.extend(sorted((ROOT / "live" / "chapters").glob("*.qmd")))

expected_tex = [
    *(ROOT / "latex" / "chapters" / f"chapter{i:02d}_{name}.tex"
      for i, name in [
          (1, "intro"), (2, "matrix"), (3, "probability"),
          (4, "regression"), (5, "rct"), (6, "causal_forest"),
          (7, "doubly_robust"), (8, "rd"), (9, "late"),
          (10, "nn"), (11, "did")]),
    ROOT / "latex" / "appendices" / "appendix_matrix.tex",
    ROOT / "latex" / "appendices" / "appendix_probability.tex",
]

for path in expected_tex:
    if not path.is_file():
        complain(path, "missing expected source")

for path in tex_files:
    text = path.read_text(encoding="utf-8")
    if re.search(r"<(font|img)\b", text, flags=re.I):
        complain(path, "contains raw HTML")
    if re.search(r"/Users/|\\Users\\", text):
        complain(path, "contains an absolute user path")
    if "mybinder.org" in text:
        complain(path, "contains an obsolete Binder link")
    if re.search(r"\$\$\s*\\begin\{(?:equation|align)", text):
        complain(path, "nests a LaTeX display environment inside $$")
    for env in ("equation", "align", "align*", "figure", "table", "rcode"):
        begins = len(re.findall(rf"\\begin\{{{re.escape(env)}\}}", text))
        ends = len(re.findall(rf"\\end\{{{re.escape(env)}\}}", text))
        if begins != ends:
            complain(path, f"unbalanced {env} environment ({begins} vs {ends})")

for path in qmd_files:
    text = path.read_text(encoding="utf-8")
    if text.count("```") % 2:
        complain(path, "unbalanced fenced code blocks")
    if re.search(r"/Users/|\\Users\\", text):
        complain(path, "contains an absolute user path")
    if "mybinder.org" in text:
        complain(path, "contains an obsolete Binder link")
    if "../../data/" in text or "../../assets/" in text:
        complain(path, "uses a path outside the live-site project")
    if 'read.csv("data/' in text and "resources:" not in text:
        complain(path, "reads course data without declaring a VFS resource")
    if re.search(r"```\s*\{?[Rr]\}?", text):
        warnings.append(f"{path.relative_to(ROOT)}: contains a non-WebR R fence")

rd_page = (ROOT / "live" / "chapters" / "chapter08_rd.qmd").read_text(encoding="utf-8")
did_page = (ROOT / "live" / "chapters" / "chapter11_did.qmd").read_text(encoding="utf-8")
for path, text in [
    (ROOT / "live" / "chapters" / "chapter08_rd.qmd", rd_page),
    (ROOT / "live" / "chapters" / "chapter11_did.qmd", did_page),
]:
    if "https://webr.r-wasm.org/v0.5.9/" not in text:
        complain(path, "missing the WebR 0.5.9 compatibility pin")

for path in [*tex_files, *qmd_files]:
    text = path.read_text(encoding="utf-8")
    for character, name in [("‑", "nonbreaking hyphen"), ("–", "en dash"), ("—", "em dash")]:
        if character in text:
            warnings.append(f"{path.relative_to(ROOT)}: contains {name}")

if errors:
    print("ERRORS")
    print("\n".join(f"- {item}" for item in errors))
if warnings:
    print("WARNINGS")
    print("\n".join(f"- {item}" for item in warnings))

if errors:
    sys.exit(1)

print(f"Checked {len(tex_files)} TeX files and {len(qmd_files)} QMD files.")
