pdflatex previewCh3 -include-directory=../ -output-directory=../
bibtex ../previewCh3 -include-directory=../
pdflatex previewCh3.tex -include-directory=../
open previewCh3.pdf
rm ../previewCh3.aux ../previewCh3.bbl ../previewCh3.blg ../previewCh3.log ../previewCh3.out ../previewCh3.pdf
