

#ifndef HYDROCALCH_INCLUDED
#define HYDROCALCH_INCLUDED

#include "global.h"

void wetfront_calc(globalpar& gp,globalvar& gv,double *v, double *deltt);

void vol_fract_calc(globalpar& gp,globalvar& gv,double *deltt);

void upbound_calc(globalvar& gv,globalpar& gp,double* deltt,std::ofstream* logPULSEfile);

#endif