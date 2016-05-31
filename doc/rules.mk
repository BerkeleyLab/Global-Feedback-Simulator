#
#Standard things

sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables

TGT_$(d)        := $(d)/user_manuals/sourcecodemap.pdf $(d)/user_manuals/structuraldocumentation.pdf $(d)/user_manuals/errata.pdf $(d)/presentations/lbnl_intro_sept13/linac_modeling_remake.pdf $(d)/presentations/lbnl_intro_sept13/linac_modeling_for_slac.pdf $(d)/publications/LLRF13/paper/system_feedback_models.pdf $(d)/publications/LLRF13/poster/poster77.pdf

$(d)/user_manuals/sourcecodemap.pdf: $(d)/diagrams/codeoutline.png
$(d)/presentations/lbnl_intro_sept13/linac_modeling_remake.pdf: $(d)/presentations/lbnl_intro_sept13/figures/ngls_feedback.pdf $(d)/presentations/lbnl_intro_sept13/figures/model_do2.pdf $(d)/presentations/lbnl_intro_sept13/figures/architecture.pdf $(d)/presentations/lbnl_intro_sept13/figures/shallowoutline.pdf $(d)/presentations/lbnl_intro_sept13/figures/TB_py_flow.pdf $(d)/presentations/lbnl_intro_sept13/figures/gen_alg.pdf $(d)/presentations/lbnl_intro_sept13/figures/question_mark2.pdf

$(d)/publications/LLRF13/paper/system_feedback_models.dvi: $(d)/publications/LLRF13/paper/system_feedback_models.tex
	cd `dirname $@` &&  latex `basename $<`
	
$(d)/publications/LLRF13/paper/system_feedback_models.ps: $(d)/publications/LLRF13/paper/system_feedback_models.dvi
	cd `dirname $@` && dvips `basename $<`

LOCAL_PS2PDF=gs -dEPSCrop -sDEVICE=pdfwrite -sOutputFile=`basename $@` -dNOPAUSE -dBATCH -dAutoRotatePages=/None -c .setpdfwrite -f `basename $<`

$(d)/publications/LLRF13/paper/system_feedback_models.pdf: $(d)/publications/LLRF13/paper/system_feedback_models.ps
	cd `dirname $@` && $(LOCAL_PS2PDF)
	
$(d)/publications/LLRF13/poster/poster77.dvi: $(d)/publications/LLRF13/poster/poster77.tex
	cd `dirname $@` &&  latex `basename $<`
	
$(d)/publications/LLRF13/poster/poster77.ps: $(d)/publications/LLRF13/poster/poster77.dvi
	cd `dirname $@` && dvips `basename $<`

$(d)/publications/LLRF13/poster/poster77.pdf: $(d)/publications/LLRF13/poster/poster77.ps
	cd `dirname $@` && $(LOCAL_PS2PDF)

$(d)/reports/physics/latex/physics_model.pdf: $(d)/reports/physics/latex/physics_model.tex
	cd `dirname $@` && pdflatex -interaction=nonstopmode `basename $<`

$(d)/reports/physics/latex/cavity_model.pdf: $(d)/reports/physics/latex/cavity_model.tex
	cd `dirname $@` && pdflatex -interaction=nonstopmode `basename $<`

LATEX_AUX_FILES := $(shell  find `pwd` -name "*.dvi" -o -name "*.aux" -o -name "*.nav" -o -name "*.toc" -o -name "*.log" -o -name "*.snm" -o -name "*.out" -type f)

CLEAN           := $(CLEAN) $(TGT_$(d)) $(d)/user_manuals/*.ps $(d)/diagrams/codeoutline.png $(LATEX_AUX_FILES) $(d)/publications/LLRF13/paper/system_feedback_models.ps $(d)/publications/LLRF13/poster/poster77.ps

# Standard things

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
