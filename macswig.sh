swig -v -python linac.i
gcc -fPIC -c linac_param.c filter.c linac_wrap.c -I/usr/include/python2.7 -I/usr/include/numpy
ld -dynamiclib linac_param.o filter.o linac_wrap.o -o _linac.so