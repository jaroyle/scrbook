(TeX-add-style-hook "EDmanuscript"
 (lambda ()
    (LaTeX-add-bibliographies
     "../AndyRefs_alphabetized")
    (LaTeX-add-labels
     "eq.encounter"
     "eq.costweighted"
     "eq.cost"
     "ecoldist.fig.raster"
     "fig.homeranges"
     "sec.mle"
     "mle.eq.cond-on-s"
     "mle.eq.intlik"
     "ecoldist.fig.raster100")
    (TeX-run-style-hooks
     "color"
     "soul"
     "verbatim"
     "hyperref"
     "amsfonts"
     "float"
     "bm"
     "natbib"
     "round"
     "colon"
     "authoryear"
     "graphicx"
     "amsmath"
     "lineno"
     "geometry"
     "8.75in}"
     "latex2e"
     "art12"
     "article"
     "12pt")))

