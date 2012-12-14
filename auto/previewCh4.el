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
     "winbugs"
     "jags"
     "secr"
     "scrbook"
     "mytilde")
    (TeX-run-style-hooks
     "fancyvrb"
     "color"
     "soul"
     "verbatim"
     "hyperref"
     "pdftex"
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "latex2e"
     "bk10"
     "book"
     "Ch4-SCR0/Chapt-SCR0")))

