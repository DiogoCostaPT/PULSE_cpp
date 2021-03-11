
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