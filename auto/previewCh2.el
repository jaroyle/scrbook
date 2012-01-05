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
    (TeX-run-style-hooks
     "verbatim"
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "latex2e"
     "bk10"
     "book"
     "Ch2/ch2")))

