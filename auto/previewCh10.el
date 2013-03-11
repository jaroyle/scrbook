(TeX-add-style-hook "previewCh10"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
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
     "pdftex"
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "bm"
     "latex2e"
     "bk10"
     "book"
     "Ch10-Design/Chapt-Design")))

