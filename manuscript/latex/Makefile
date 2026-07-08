TEX  = AllelePerturb_manuscript
BIB  = references

.PHONY: all clean

all: $(TEX).pdf

$(TEX).pdf: $(TEX).tex $(BIB).bib
	pdflatex -interaction=nonstopmode $(TEX)
	bibtex $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

clean:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out \
	      $(TEX).toc $(TEX).lof $(TEX).lot $(TEX).fls $(TEX).fdb_latexmk \
	      $(TEX).synctex.gz
