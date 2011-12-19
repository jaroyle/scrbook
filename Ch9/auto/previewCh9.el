(TeX-add-style-hook "previewCh9"
 (lambda ()
    (LaTeX-add-bibliographies
     "../AndyRefs_alphabetized")
    (TeX-run-style-hooks
     "natbib"
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "float"
     "latex2e"
     "bk10"
     "book"
     "Ch9")))

