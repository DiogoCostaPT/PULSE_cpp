

#ifndef TOOLBOXH_INCLUDED
#define TOOLBOXH_INCLUDED

#include<armadillo>
#include <vector> 
#include<iostream>
#include <dirent.h>

#include "outwrite.h"

int findLastStep(const char *path);

void checkmesh2(int* H_local,int* L_local,int* h_layer,int* l_layer,int* nh,int* nl,std::ofstream* logPULSEfile);

#endif