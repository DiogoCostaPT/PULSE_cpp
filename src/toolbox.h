

#ifndef TOOLBOXH_INCLUDED
#define TOOLBOXH_INCLUDED

#include<armadillo>
#include <vector> 
#include<iostream>
#include <dirent.h>

#include "outwrite.h"

std::string SplitFilename (const std::string& str);

int findLastStep(const char *path);

void checkmesh2(double* H_local,double* L_local,double* h_layer,double* l_layer,
                int* nh,int* nl,std::ofstream* logPULSEfile,std::string* results_flname);

bool read_matrixes_ext(globalpar& gp,globalvar& gv,
            std::string* time_file, std::string* v_ice_file,std::string* v_liquid_file,std::string* v_ice2liq_1_file,
            std::string* v_ice2liq_2_file, std::string*fluxQ_file,std::ofstream* logPULSEfile);

std::string removeSpaces(std::string str) ;

#endif