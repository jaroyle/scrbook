(TeX-add-style-hook "previewCh2"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (LaTeX-add-labels
     "chapt.intro"
     "chapt.closed"
     "chapt.scr0"
     "chapt.poisson"
     "chapt.mle"
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
     "Ch2/Ch2")))

