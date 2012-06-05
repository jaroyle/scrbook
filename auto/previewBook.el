(TeX-add-style-hook "previewBook"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (TeX-add-symbols
     "R"
     "winbugs"
     "bugs"
     "jags"
     "secr"
     "scrbook")
    (TeX-run-style-hooks
     "url"
     "bm"
     "makeidx"
     "fancyvrb"
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
     "elsst-book"
     "latex2e"
     "bk10"
     "book"
     "Ch2/Ch2"
     "Ch3/Ch3"
     "Ch4/Chapter4"
     "Ch9/Ch9"
     "Ch7/Ch7")))

