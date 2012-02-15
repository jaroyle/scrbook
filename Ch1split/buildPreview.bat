pdflatex previewCh1split -include-directory=../ -output-directory=../
bibtex ../previewCh1split -include-directory=../
pdflatex previewCh1split.tex -include-directory=../
open previewCh1split.pdf
rm ../previewCh1split.aux ../previewCh1split.bbl ../previewCh1split.blg ../previewCh1split.log ../previewCh1split.out ../previewCh1split.pdf
