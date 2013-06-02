pdflatex previewCh2 -include-directory=../ -output-directory=../
bibtex ../previewCh2 -include-directory=../
pdflatex previewCh2.tex -include-directory=../
open previewCh2.pdf
rm ../previewCh2.aux ../previewCh2.bbl ../previewCh2.blg ../previewCh2.log ../previewCh2.out ../previewCh2.pdf
