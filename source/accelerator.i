%module accelerator

%{
#include "filter.h"
#include "cavity.h"
#include "rf_station.h"
%}

%include "complex.i"
%include "carrays.i"
%array_class(int, intArray);
%array_class(double complex, complexdouble_Array);
%array_class(double, double_Array);

%include "cpointer.i"
%pointer_class(int, intp);
%pointer_class(double complex, compp);

%include "filter.h"
%include "cavity.h"
%include "rf_station.h"