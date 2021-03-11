
// Copyright 2021: Diogo Costa

// This program, PULSE_cpp, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
            std::string* time_file, std::string* v_ice_file,std::string* v_liquid_file,
            std::string* v_ice2liq_1_file,std::string* v_ice2liq_2_file, std::string*fluxQ_file, 
            std::string* prec_c_ext_file,std::ofstream* logPULSEfile);

std::string removeSpaces(std::string str) ;

#endif