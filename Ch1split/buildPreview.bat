pdflatex previewCh1 -include-directory=../ -output-directory=../
bibtex ../previewCh1 -include-directory=../
pdflatex previewCh1.tex -include-directory=../
open previewCh1.pdf
rm ../previewCh1.aux ../previewCh1.bbl ../previewCh1.blg ../previewCh1.log ../previewCh1.out ../previewCh1.pdf
