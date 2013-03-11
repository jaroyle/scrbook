(TeX-add-style-hook "previewCh14"
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
     "fancyvrb"
     "latex2e"
     "bk10"
     "book"
     "Ch14-Multisession/Chapt-Multisession")))

