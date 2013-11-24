#Build flags for all targets

GCC_FLAGS = -Wstrict-prototypes -Wpointer-arith -Wcast-align -Wcast-qual \
    -Wshadow -Waggregate-return -Wmissing-prototypes -Wnested-externs \
    -Wall -W -Wno-unused -Winline -Wwrite-strings -Wundef -pedantic
CF_ALL = -Wall -O2 -W -fPIC -g -std=c99 -D_GNU_SOURCE $(GCC_FLAGS) ${CFLAGS_$@}
LF_ALL = -lm
LL_ALL =

OCTAVE = octave
AWK= awk
LYX = lyx
DOT = dot
PDFLATEX = pdflatex
DVIPS = dvips
LATEX = latex

#Build tools

CC = ./build/ccd-gcc
INST = ./build/install
COMP = $(CC) $(CF_ALL) $(CF_TGT) -o $@ -c $< ${CFLAGS_$@}
LINK = $(CC) $(LF_ALL) $(LF_TGT) -o $@ $^ $(LL_TGT) $(LL_ALL)
COMPLINK = $(CC) $(CF_ALL) $(CF_TGT) $(LF_ALL) $(LF_TGT) -o $@ $< $(LL_TGT) $(LL_ALL)
ARCH = ar rcs $@ $^
LYX2PS = $(LYX) --export ps $<
DOT2PNG = $(DOT) -Tpng -o $@ $<

OCTAVE_SILENT=$(OCTAVE) -q $<
PS2PDF=gs -dEPSCrop -sDEVICE=pdfwrite -sOutputFile=$@ -dNOPAUSE -dBATCH -dAutoRotatePages=/None -c .setpdfwrite -f $<
TEX2PDF = cd `dirname $@` && $(PDFLATEX) `basename $<`
DVI2PS = cd `dirname $@` && $(DVIPS) `basename $<`

#Standard parts

include rules.mk
