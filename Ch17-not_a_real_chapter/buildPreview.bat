pdflatex previewCh17 -include-directory=../ -output-directory=../
bibtex ../previewCh17 -include-directory=../
pdflatex previewCh17.tex -include-directory=../
open previewCh17.pdf
rm ../previewCh17.aux ../previewCh17.bbl ../previewCh17.blg ../previewCh17.log ../previewCh17.out ../previewCh17.pdf
