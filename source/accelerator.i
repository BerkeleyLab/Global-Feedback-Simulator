%module accelerator

%{
#define SWIG_FILE_WITH_INIT
#include "filter.h"
#include "cavity.h"
#include "rf_station.h"
#include "cryomodule.h"
#include "linac.h"
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
%include "cryomodule.h"
%include "linac.h"
