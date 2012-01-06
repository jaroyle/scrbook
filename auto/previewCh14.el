(TeX-add-style-hook "previewCh14"
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
    (TeX-run-style-hooks
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
     "Ch14/Ch14")))

