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

#ifndef OUTWRITEH_INCLUDED
#define OUTWRITEH_INCLUDED

#include<iostream>
#include<armadillo>

#include "global.h"

void print_screen_log(std::ofstream* logPULSEfile,std::string* msg);

bool print_results(globalvar& gv,globalpar& gp, int print_tag, unsigned int print_step, 
        std::chrono::duration<double> elapsed_seconds,std::string* results_flname);

#endif