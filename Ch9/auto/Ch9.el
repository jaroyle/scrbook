(TeX-add-style-hook "Ch9"
 (lambda ()
    (LaTeX-add-labels
     "eq:pdf-ipp"
     "fig:elevMap"
     "tab:simIPP")
    (TeX-run-style-hooks
     "lineno"
     "graphicx"
     "amsfonts"
     "amsmath"
     "latex2e"
     "bk10"
     "book")))

