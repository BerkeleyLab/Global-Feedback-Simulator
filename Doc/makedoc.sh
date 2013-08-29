#!/bin/sh

dot -Tpng -o diagrams/codeoutline.png diagrams/codeoutline.dot 

lyx --export ps structuraldocumentation.lyx
lyx --export ps sourcecodemap.lyx
lyx --export ps errata.lyx
