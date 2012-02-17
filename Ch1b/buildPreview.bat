pdflatex previewCh1b -include-directory=../ -output-directory=../
bibtex ../previewCh1b -include-directory=../
pdflatex previewCh1b.tex -include-directory=../
open previewCh1b.pdf
rm ../previewCh1b.aux ../previewCh1b.bbl ../previewCh1b.blg ../previewCh1b.log ../previewCh1b.out ../previewCh1b.pdf
