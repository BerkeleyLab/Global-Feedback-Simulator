#
#Standard things

sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables

TGT_$(d)        := $(d)/reports/physics/latex/physics_model.pdf

$(d)/reports/physics/latex/physics_model.pdf: $(d)/reports/physics/latex/physics_model.tex
	cd `dirname $@` && pdflatex -interaction=nonstopmode `basename $<`

$(d)/reports/physics/latex/cavity_model.pdf: $(d)/reports/physics/latex/cavity_model.tex
	cd `dirname $@` && pdflatex -interaction=nonstopmode `basename $<`

LATEX_AUX_FILES := $(shell  find `pwd` -name "*.dvi" -o -name "*.aux" -o -name "*.nav" -o -name "*.toc" -o -name "*.log" -o -name "*.snm" -o -name "*.out" -type f)

CLEAN           := $(CLEAN) $(TGT_$(d)) $(LATEX_AUX_FILES)

# Standard things

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
