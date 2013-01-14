(TeX-add-style-hook "previewCh11"
 (lambda ()
    (LaTeX-add-bibliographies
     "../AndyRefs_alphabetized")
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
     "Ch11")))

