(TeX-add-style-hook "previewCh4"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (LaTeX-add-labels
     "chapt.intro"
     "chapt.glms"
     "chapt.closed"
     "chapt.mle"
     "chapt.mcmc"
     "chapt.ipp")
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
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "latex2e"
     "bk10"
     "book"
     "Chapter4")))

