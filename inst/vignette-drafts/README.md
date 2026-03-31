# Draft vignettes (not built on pkgdown)

These `.Rmd` files are **not** under `vignettes/`, so they are **not** rendered when building the pkgdown site or when `R CMD check` builds vignettes. That keeps CI and the documentation site stable while still shipping the sources in the repository.

**FSHD (Scissor-style)** and **Alzheimer’s (custom reference)** drafts live here with any small sidecar files (e.g. CyteTypeR annotation CSV).

Render locally from the package root, for example:

```r
rmarkdown::render("inst/vignette-drafts/fshd-scissor-custom-reference.Rmd")
rmarkdown::render("inst/vignette-drafts/alzheimers-custom-reference-draft.Rmd")
```

To publish them on the site again later, move the `.Rmd` back into `vignettes/` and re-add entries to `_pkgdown.yml` once downloads, data, and runtime are reliable in CI.
