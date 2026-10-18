# Build the GLoBES.jl manual: `make` computes the figure data with GLoBES.jl and
# typesets the manual with LuaLaTeX.

JULIA ?= julia --project=../..
LATEX ?= lualatex -interaction=nonstopmode -halt-on-error

all: figures pdf

figures: figures/data/tour-output.txt figures/data/prob_nu.dat

figures/data/prob_nu.dat: figures/make_figures.jl
	$(JULIA) figures/make_figures.jl

figures/data/tour-output.txt: ../../examples/tour.jl
	mkdir -p figures/data
	$(JULIA) ../../examples/tour.jl > $@

pdf:
	$(LATEX) globes-jl-manual.tex
	$(LATEX) globes-jl-manual.tex

clean:
	rm -f *.aux *.log *.out *.toc chapters/*.aux

.PHONY: all figures pdf clean
