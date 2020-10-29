
#ifndef INITIATEH_INCLUDED
#define INITIATEH_INCLUDED

#include "global.h"
#include "toolbox.h"
#include "outwrite.h"

bool initiate(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile,
    std::string* results_flname);

#endif