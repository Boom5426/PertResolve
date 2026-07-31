TEX  = AllelePerturb_manuscript
SI   = AllelePerturb_SI
BIB  = references

# The figures are real prerequisites: without them `make` reports success while
# leaving the previous figure PDFs embedded in the output, so a rebuilt panel
# silently does not reach the manuscript.
FIGS = $(wildcard figures/*.pdf)

.PHONY: all si clean

all: $(TEX).pdf $(SI).pdf

$(TEX).pdf: $(TEX).tex $(BIB).bib $(FIGS)
	pdflatex -interaction=nonstopmode $(TEX)
	bibtex $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

# Supplementary Information (no bibliography; twice resolves table labels)
si: $(SI).pdf

$(SI).pdf: $(SI).tex $(FIGS)
	pdflatex -interaction=nonstopmode $(SI)
	pdflatex -interaction=nonstopmode $(SI)

clean:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out \
	      $(TEX).toc $(TEX).lof $(TEX).lot $(TEX).fls $(TEX).fdb_latexmk \
	      $(TEX).synctex.gz \
	      $(SI).aux $(SI).log $(SI).out $(SI).toc $(SI).lot $(SI).lof \
	      $(SI).fls $(SI).fdb_latexmk $(SI).synctex.gz
