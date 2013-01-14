pdflatex previewCh11 -include-directory=../ -output-directory=../
bibtex ../previewCh11 -include-directory=../
pdflatex previewCh11.tex -include-directory=../
open previewCh11.pdf
rm ../previewCh11.aux ../previewCh11.bbl ../previewCh11.blg ../previewCh11.log ../previewCh11.out ../previewCh11.pdf
