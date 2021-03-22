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

#include "hydroinput.h"

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpMeteo(globalvar& gv,double *tcum)
{
    double meteo_t = 0.0f,meteo_t2 = 0.0f, // time
            temp_calc_t = 0.0f,
            rainfall_calc_t = 0.0f,
            snowfall_calc_i = 0.0f, // melt rate
            precipt_c_calc_t = 0.0f; // concentration of added snow
    unsigned int nprec = 0;;
    
    nprec = int((*gv.meteoall_int).n_rows) - 1;
    
    // identification of the time step and linear interpolation for time t
    meteo_t2 = 0.0f;
    for(int a=0;a<nprec;a++){
        meteo_t = (*gv.meteoall_int).at(a,0);
        meteo_t2 = (*gv.meteoall_int).at(a+1,0);
        
        if ((*tcum) > meteo_t && (*tcum) <= meteo_t2) {
            if ((meteo_t2 - meteo_t) > 0){
                temp_calc_t = (*gv.meteoall_int).at(a,1);
                rainfall_calc_t = (*gv.meteoall_int).at(a,2);
                snowfall_calc_i = (*gv.meteoall_int).at(a,3);
                precipt_c_calc_t = (*gv.meteoall_int).at(a,4);
                
                gv.tempert_t = temp_calc_t;
                gv.rainfall_t = rainfall_calc_t / (meteo_t2 - meteo_t); // m/deltatime -> m/sec
                gv.snowfall_t = snowfall_calc_i / (meteo_t2 - meteo_t); // m/deltatime -> m/sec
                gv.precip_c_t = precipt_c_calc_t;        
             }else
            {
                gv.snowfall_t = 0.0f;
            }
        break;
        }
    
    }
}

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpQmelt(globalvar& gv,double *tcum)
{

    double qmlt_time1 = 0.0f,qmlt_time2 = 0.0f, // time
            qmlt_t = 0.0f; // melt rate
    unsigned int nqmlt = 0;;
    
    nqmlt = int((*gv.meteoall_int).n_rows) - 1;
    
    // identification of the time step and linear interpolation for time t
    qmlt_time2 = 0.0f;
    for(int a=0;a<nqmlt;a++){
        qmlt_time1 = (*gv.qcmel_int).at(a,0);
        qmlt_time2 = (*gv.qcmel_int).at(a+1,0);
        qmlt_t = (*gv.qcmel_int).at(a,1);

        if ((*tcum) > qmlt_time1 && (*tcum) <= qmlt_time2) {
            if ((qmlt_time2 - qmlt_time1) > 0){
                gv.qmelt_t = qmlt_t / (qmlt_time2 - qmlt_time1); // m/deltatime -> m/sec      
             }else
            {
                gv.qmelt_t = 0.0f;
            }
        break;
        }

    }
}