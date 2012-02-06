pdflatex previewCh6 -include-directory=../ -output-directory=../
bibtex ../previewCh6 -include-directory=../
pdflatex previewCh6.tex -include-directory=../
open previewCh6.pdf
rm ../previewCh6.aux ../previewCh6.bbl ../previewCh6.blg ../previewCh6.log ../previewCh6.out ../previewCh6.pdf
