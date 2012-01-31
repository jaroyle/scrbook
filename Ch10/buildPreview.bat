pdflatex previewCh10 -include-directory=../ -output-directory=../
bibtex ../previewCh10 -include-directory=../
pdflatex previewCh10.tex -include-directory=../
open previewCh10.pdf
rm ../previewCh10.aux ../previewCh10.bbl ../previewCh10.blg ../previewCh10.log ../previewCh10.out ../previewCh10.pdf
