

#ifndef HYDROCALCH_INCLUDED
#define HYDROCALCH_INCLUDED

#include "global.h"

void vol_fract_calc(globalpar& gp,globalvar& gv,double *v,double *deltt);

void wetfront_calc(globalpar& gp,globalvar& gv,double *v, double *deltt);

void watermass_calc_internal(globalvar& gv,globalpar& gp,double* deltt,double *v,
    std::ofstream* logPULSEfile);

bool watermass_calc_external(globalvar& gv,globalpar& gp,double* deltt,
        std::ofstream* logPULSEfile, int t);

#endif