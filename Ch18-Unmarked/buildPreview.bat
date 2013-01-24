pdflatex previewCh14 -include-directory=../ -output-directory=../
bibtex ../previewCh14 -include-directory=../
pdflatex previewCh14.tex -include-directory=../
open previewCh14.pdf
rm ../previewCh14.aux ../previewCh14.bbl ../previewCh14.blg ../previewCh14.log ../previewCh14.out ../previewCh14.pdf
