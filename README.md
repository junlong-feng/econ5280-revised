# ECON 5280 revised lecture notes

This project contains two coordinated versions of the course:

1. `latex/`: the canonical Palatino lecture-note book and printable PDF.
2. `live/`: Quarto Live pages whose R cells run in the browser through WebR.

The original ZIP Markdown files are retained unchanged in `source_md/`. The
revised sources correct mathematical and typographical errors, place optional
matrix/probability material in appendices, and add a modern DID/event-study
chapter.

## Build the PDF

With TeX Live and `latexmk` installed:

```bash
make pdf
```

The finished file is written to
`output/pdf/ECON5280_Revised_Lecture_Notes.pdf`.

## Build or preview the live site

The Quarto Live extension is vendored in `live/_extensions/`, so the classroom
computer does not need R and the build does not need to download an extension.
With Quarto installed, preview the site from a terminal:

```bash
cd live
quarto preview
```

For a noninteractive build:

```bash
make live
```

The generated site is in `live/_site/` and can be hosted on GitHub Pages,
Quarto Pub, Netlify, or any other static web host.

Do not open the generated HTML through a `file://` URL. Preview it through
Quarto or publish `_site/` to a web server so that the browser can load the
WebAssembly runtime, packages, and course CSV files correctly.

## Classroom workflow

- Open the relevant live chapter shortly before the code demonstration.
- Run the first small cell before class to warm the WebR runtime and package
  cache in that browser.
- Use the browser cells for short examples.
- Use the included Codespaces configuration for full bootstrap inference or
  computationally intensive forest fits.
- Stop the Codespace after the demonstration; Quarto Live itself has no cloud
  session to stop.

The included GitHub Pages workflow builds and publishes the static live site
without executing R. See `runtime/README.md` for version pins, package limits,
and the native-R fallback.

## What changed

- The canonical printable notes now use Palatino text and mathematics.
- Optional matrix-algebra and probability details are retained in Appendices A
  and B rather than removed.
- Every chapter has a Quarto Live companion page.
- Chapter 11 adds modern difference-in-differences and event-study methods.
- The missing Chapter 8 figure was reconstructed reproducibly from `house.csv`.

See `REVISION_NOTES.md` for the chapter-by-chapter correction log.
