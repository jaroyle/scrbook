pdflatex previewCh5 -include-directory=../ -output-directory=../
bibtex ../previewCh5 -include-directory=../
pdflatex previewCh5.tex -include-directory=../
open previewCh5.pdf
rm ../previewCh5.aux ../previewCh5.bbl ../previewCh5.blg ../previewCh5.log ../previewCh5.out ../previewCh5.pdf
