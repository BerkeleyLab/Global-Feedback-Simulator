# Standard things

sp := $(sp).x
dirstack_$(sp) := $(d)
d := $(dir)

# Local variables

#TGT_$(d) :=
#TGT_$(d)_TB :=

#CHK_$(d) :=

# Local variables

CLEAN := $(CLEAN) $(TGT_$(d))

# Standard things

-include $(DEPS_$(d))

d  := $(dirstack_$(sp))
sp := $(basename $(sp))
