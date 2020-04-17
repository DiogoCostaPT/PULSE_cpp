
#ifndef OUTWRITEH_INCLUDED
#define OUTWRITEH_INCLUDED

#include<iostream>
#include<armadillo>

#include "global.h"

void print_screen_log(std::ofstream* logPULSEfile,std::string* msg);

bool print_results(globalvar& gv,globalpar& gp, int print_tag, unsigned int print_step, 
        std::chrono::duration<double> elapsed_seconds,std::string* results_flname);

#endif