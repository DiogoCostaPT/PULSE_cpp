
#ifndef PULSEH__INCLUDED
#define PULSEH__INCLUDED

#include <armadillo>

#include "global.h"
#include "toolbox.h"
#include "hydroinput.h"
#include "hydrocalc.h"
#include "crank_nicholson_solve.h"

void pulsemodel(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile);

#endif