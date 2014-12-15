#
#Standard things

sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables

TGT_$(d)        := 

$(d)/linac.py: $(d)/linac.i
	swig -v -python $^

$(d)/linac_wrap.c: $(d)/linac.i
	swig -v -python $^

CFLAGS_$(d)/linac_param.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/doublecompress.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/filter.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/step_llrf.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/state_space_top.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/beam_based_feedback.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/dynamic_noise.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/cavity.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/linac_wrap.o := -I/usr/include/python2.7 -I/usr/include/numpy

$(d)/_linac.so: $(d)/linac_param.o $(d)/doublecompress.o $(d)/filter.o $(d)/cavity.o $(d)/step_llrf.o $(d)/state_space_top.o $(d)/beam_based_feedback.o $(d)/dynamic_noise.o $(d)/linac_wrap.o
	ld -shared $^ -o $@

CLEAN           := $(CLEAN) $(TGT_$(d)) $(d)/linac.py $(d)/linac_wrap.c $(d)/_linac.so $(d)/filter.py $(d)/*.o $(d)/*.o.d $(d)/*.pyc

# Standard things

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
