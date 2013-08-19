#!/usr/bin/python

# for now... lets use a makefile later
import os
os.system('gcc -fPIC -c filter.c -o filter.so')

from ctypes import *

lib_filter = CDLL('filter.so')

class Filter(Structure):
    _fields_ = [("alloc_order",c_int),
                ("alloc_coeffs",c_int),
                ("n_coeffs",c_int),
                ("order",c_int),
                ("modes",POINTER(c_int)),
                ("coeff_start",POINTER(coeff_start)),
                ("coeffs",POINTER(c_motherfuckernocomplex)]
