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

CFLAGS_$(d)/filter.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/rf_station.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/cryomodule.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/linac.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/doublecompress.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/simulation_top.o := -I/usr/include/python2.7 -I/usr/include/numpy
# CFLAGS_$(d)/beam_based_feedback.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/noise.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/cavity.o := -I/usr/include/python2.7 -I/usr/include/numpy
CFLAGS_$(d)/accelerator_wrap.o := -I/usr/include/python2.7 -I/usr/include/numpy

$(d)/accelerator_wrap.o: CF_ALL := $(filter-out -Wcast-qual -Wshadow -Wmissing-prototypes -Wstrict-prototypes,$(CF_ALL))
$(d)/_accelerator.so: $(d)/filter.o $(d)/cavity.o $(d)/rf_station.o $(d)/cryomodule.o $(d)/linac.o $(d)/doublecompress.o $(d)/noise.o $(d)/simulation_top.o $(d)/accelerator_wrap.o
	$(CC) -shared $^ -o $@
	# Use this rule if you're running under Mac OSX
	# $(CC) -lpython -dynamclib $^ -o $@

export UNIT_TEST_FILES := $(d)/unit_tests_all.py $(d)/cavity_test.py $(d)/rf_station_test.py $(d)/cryomodule_test.py $(d)/doublecompress_test.py $(d)/simulation_test.py

CLEAN           := $(CLEAN) $(TGT_$(d)) $(d)/accelerator.py $(d)/accelerator_wrap.c $(d)/_accelerator.so $(d)/*.o $(d)/*.o.d $(d)/*.pyc

# Standard things

d               := $(dirstack_$(sp))
sp              := $(basename $(sp))
