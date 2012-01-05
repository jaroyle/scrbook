pdflatex previewCh4 -include-directory=../ -output-directory=../
bibtex ../previewCh4 -include-directory=../
pdflatex previewCh4.tex -include-directory=../
open previewCh4.pdf
rm ../previewCh4.aux ../previewCh4.bbl ../previewCh4.blg ../previewCh4.log ../previewCh4.out ../previewCh4.pdf
