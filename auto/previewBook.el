(TeX-add-style-hook "previewBook"
 (lambda ()
    (LaTeX-add-bibliographies
     "AndyRefs_alphabetized")
    (TeX-add-symbols
     "mytilde"
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
     "rotating"
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
     "elsst-book"
     "alltt"
     "latex2e"
     "bk10"
     "book"
     "Ch1/Chapt1"
     "Ch2/Chapt2"
     "Ch3-Bayes/Chapt-Bayes"
     "Ch4-Closed/Chapt-Closed"
     "Ch5-SCR0/Chapt-SCR0"
     "Ch6-MLE/Chapt-MLE"
     "Ch7-Covariates/Chapt-Covariates"
     "Ch8-GoF/Chapt-GoF"
     "Ch9-PoisMn/Chapt-PoisMn"
     "Ch10-Design/Chapt-Design"
     "Ch11-Statespace/Chapt-Statespace"
     "Ch12-EcolDist/Chapt-EcolDist"
     "Ch13-RSF/Chapt-RSF"
     "Ch14-Multisession/Chapt-Multisession"
     "Ch15-searchencounter/Chapt-searchencounter"
     "Ch16-Open/Chapt-Open"
     "Ch17-MCMC/Chapt-MCMC"
     "Ch18-Unmarked/Chapt-Unmarked"
     "Ch19-PartialID/Chapt-PartialID"
     "Ch20-Last/Chapt-Last"
     "AppendixI/AppendixI")))

