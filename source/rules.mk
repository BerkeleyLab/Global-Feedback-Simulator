#
#Standard things

sp              := $(sp).x
dirstack_$(sp)  := $(d)
d               := $(dir)

# Local variables

TGT_$(d)        := 

$(d)/accelerator.py: $(d)/accelerator.i
	swig -v -python $^

$(d)/accelerator_wrap.c: $(d)/accelerator.i
	swig -v -python $^

# CFLAGS_$(d)/linac_param.o := -I/usr/include/python2.7 -I/usr/include/numpy
# CFLAGS_$(d)/doublecompress.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/filter.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/rf_station.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/cryomodule.o := -I/usr/include/python2.7 -I/usr/include/numpy
# CFLAGS_$(d)/state_space_top.o := -I/usr/include/python2.7 -I/usr/include/numpy
# CFLAGS_$(d)/beam_based_feedback.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/dynamic_noise.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/cavity.o := -I/usr/include/python2.7 -I/usr/include/numpy
# CFLAGS_$(d)/linac_wrap.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/accelerator_wrap.o := -I/usr/include/python2.7 -I/usr/include/numpy

# $(d)/_linac.so: $(d)/linac_param.o $(d)/doublecompress.o $(d)/filter.o $(d)/cavity.o $(d)/step_llrf.o $(d)/state_space_top.o $(d)/beam_based_feedback.o $(d)/dynamic_noise.o $(d)/linac_wrap.o
	# ld -shared $^ -o $@

$(d)/accelerator_wrap.o: CF_ALL := $(filter-out -Wcast-qual -Wshadow -Wmissing-prototypes,$(CF_ALL))
$(d)/_accelerator.so: $(d)/filter.o $(d)/cavity.o $(d)/rf_station.o $(d)/cryomodule.o $(d)/accelerator_wrap.o
	$(CC) -shared $^ -o $@

export UNIT_TEST_FILES := $(d)/unit_tests_all.py $(d)/cavity_test.py $(d)/rf_station_test.py $(d)/cryomodule_test.py

CLEAN           := $(CLEAN) $(TGT_$(d)) $(d)/accelerator.py $(d)/accelerator_wrap.c $(d)/_accelerator.so $(d)/*.o $(d)/*.o.d $(d)/*.pyc

# Standard things

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
