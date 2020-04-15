

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
                int* nh,int* nl,std::ofstream* logPULSEfile);

#endif