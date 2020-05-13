
#ifndef HYDROINPUTH_INCLUDED
#define HYDROINPUTH_INCLUDED


#include<armadillo>
#include "global.h"

void findInterpMeteo(globalvar& gv,double *tcum);

void findInterpQmelt(globalvar& gv,double *tcum);

#endif