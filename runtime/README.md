# Runtime design

The primary classroom runtime is Quarto Live with WebR. Each chapter declares
only the packages it needs. Most chapters use the current WebR engine. The RD
and DID pages pin WebR 0.5.9 (R 4.5) because that binary repository contains
`rdd` and the complete `did`/`fastglm`/`fixest` dependency chain. Keeping the
pin in those two pages avoids a package-resolution failure if the current WebR
repository has a different package set.

WebR can load only packages compiled for WebAssembly. It cannot compile an
ordinary CRAN or GitHub source package in the browser. Browser examples should
therefore remain single-threaded and modest in size.

The `.devcontainer/devcontainer.json` file provides a two-core native R 4.5
Codespace as a fallback. Use it for full `grf` fits, the full 1,000-iteration
multiplier bootstrap, or any package whose WebAssembly build is unavailable.

## Recommended classroom sequence

1. Publish `live/_site/` to a static HTTPS host such as GitHub Pages.
2. Open the day's chapter before class and run its first small cell once.
3. Keep that tab open during the lecture; later cells reuse the same in-browser
   R session.
4. Use the Codespace only for the explicitly marked full-size computations.

The live examples deliberately use smaller forests, samples, and bootstrap
counts. Those changes affect runtime only; the LaTeX notes state the settings
recommended for substantive inference.

## Rebuilding

- Quarto CLI: 1.10.18 (also pinned in the dev container and Pages workflow)
- Native fallback: Rocker R 4.5 dev container
- Browser runtime: current WebR except for the RD and DID pins noted above

The Quarto Live extension is included under `live/_extensions/r-wasm/live/`,
along with its license, so a fresh extension download is not required.
