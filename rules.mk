#Standard things

.SUFFIXES:
.SUFFIXES: .c .o

all: targets
check: check_all

#Subdirectories in random order,

include dir_list.mk
dir := $(BUILD_DIR)
include     $(dir)/rules.mk
dir := $(DOC_DIR)
include     $(dir)/rules.mk
dir := $(OUTPUT_DIR)
include     $(dir)/rules.mk
dir := $(PLOTTING_DIR)
include     $(dir)/rules.mk
dir := $(SOURCE_DIR)
include     $(dir)/rules.mk

#General directory-independent rules

%.o:    %.c
		$(COMP)

%:      %.o
		$(LINK)

%:      %.c
		$(COMPLINK)

%.a:    %.o
		$(ARCH)

%_check:    %_tb build/testcode.awk
		$(VERILOG_CHECK)

%:      %.m
		$(OCTAVE_SILENT)

%.pdf:  %.eps
		$(PS2PDF)

%.pdf:  %.ps
		$(PS2PDF)

%.ps:	%.lyx
		$(LYX2PS)

%.png:	%.dot
		$(DOT2PNG)

%.pdf:	%.tex
		$(TEX2PDF)

CLEAN := $(CLEAN) *.log physics_model.pdf cavity_model.pdf

# The variables TGT_*, CLEAN and CMD_INST* may be added to by the Makefile
# fragments in the various subdirectories

ePHONY:     targets
targets:    $(TGT_$(DOC_DIR)) $(TGT_$(SOURCE_DIR)) $(TGT_$(UNIT_TESTS_DIR))

.PHONY:     check_all
check_all:      $(CHK_ALL_$(SOURCE_DIR)) $(CHK_ALL_$(UNIT_TESTS_DIR))

# The "find" commands below check that the source code satisifies:
#  no hidden files
#  no filenames with spaces
#  files don't contain trailing spaces or tabs (.eps files excepted)
.PHONY:     clean
clean:
		rm -f $(CLEAN)
#		rm -rf $(CLEAN_DIRS)
#		! find . -path ./.git -prune -o -name ".*" -and -not -name "." -and -not -name ".gitignore" -print | grep .
#		! find . -path ./.git -prune -o -print | grep " "
#		! find . -path ./.git -prune -o -type f -and -not -name "*.eps" -and -not -name "*.ps" -and -not -name "*.gold" -print | LC_ALL=C xargs grep -n -E "`printf \"\t$$| $$| \t|[^[:alnum:][:punct:] \t]\"`"

# Prevent make from removing any build targets, including intermediate ones

.SECONDARY: $(CLEAN)

unit_tests.log: $(SOURCE_DIR)/_accelerator.so $(UNIT_TEST_FILES)
	python $(SOURCE_DIR)/unit_tests_all.py > $@

physics_model.pdf: $(DOC_DIR)/reports/physics/latex/physics_model.pdf
	mv $< $@

cavity_model.pdf: $(DOC_DIR)/reports/physics/latex/cavity_model.pdf
	mv $< $@
