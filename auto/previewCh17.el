(TeX-add-style-hook "previewCh17"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (TeX-run-style-hooks
     "color"
     "soul"
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
     "Ch17/Ch17")))

