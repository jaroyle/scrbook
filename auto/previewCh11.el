(TeX-add-style-hook "previewCh11"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (LaTeX-add-labels
     "chapt.intro"
     "chapt.glms"
     "chapt.closed"
     "chapt.scr0"
     "chapt.poisson"
     "chapt.mcmc"
     "chapt.gof"
     "chapt.covariates"
     "chapt.ipp"
     "chapt.open")
    (TeX-add-symbols
     "R"
     "bugs"
     "jags"
     "secr"
     "scrbook")
    (TeX-run-style-hooks
     "color"
     "soul"
     "verbatim"
     "hyperref"
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "latex2e"
     "bk10"
     "book"
     "Ch11/Ch11")))

