
// Copyright 2017-2019: Diogo Costa

// This program, PULSE, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) an_col later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include<iostream>
#include<fstream>
#include<math.h>
#include<armadillo>
#include<string>
#include<memory> 
#include <chrono>
#include <ctime>  

#include <vector> 
#include <dirent.h>
#include <sys/types.h>

#include "global.h"
#include "outwrite.h"
#include "toolbox.h"
#include "hydroinput.h"
#include "readin.h"
#include "hydrocalc.h"
#include "crank_nicholson_solve.h"
#include "pulse.h"
#include "initiate.h"


int main(int argc, char* argv[]) 
{   
    
    double H_local,L_local,h_layer,l_layer,vfrac_air_frshsnow;
    int nl,nh;
    int n_qcmelt,n_meteoall;
    std::string sim_purp,qcmelt_file,meteo_file,msg;
    std::ofstream logPULSEfile ("log.pulse");
    
    std::string modset_flname (argv[1]);
    std::string results_flname (argv[2]);

    // Assign global parameters
    globalpar gp; 
    
    //try{
    
        msg = "......................................... \n PULSE: multi-phase multi-layer snowpack chemistry model \n......................................... \n";
        print_screen_log(&logPULSEfile,&msg);   

        // read simulation setup
        read_simset(gp,modset_flname,&sim_purp,
            &h_layer,&l_layer,&qcmelt_file,
            &meteo_file,&logPULSEfile,&n_qcmelt,&n_meteoall,&vfrac_air_frshsnow);  

        // Check if input data is consistent
        if (n_qcmelt!=n_meteoall){
            msg = "Input timeseries need to have the same size: " + (qcmelt_file)
                + " and " + (meteo_file);   
            print_screen_log(&logPULSEfile,&msg); 
            abort();
        }

        // create mesh
        //checkmesh(&H_local,&L_local,&h_layer,&l_layer,&nh,&nl,&logPULSEfile);
        checkmesh2(&H_local,&L_local,&h_layer,&l_layer,&nh,&nl,&logPULSEfile,&results_flname);

        // Asign global variables (heap)
        globalvar gv(nh,nl,n_qcmelt,n_meteoall); 
        gv.snowH = H_local;
        gv.snowL = L_local;
        gv.snowh = h_layer;
        gv.snowl = l_layer;
        gv.vfrac_air_frshsnow = vfrac_air_frshsnow;

        // read snowmelt input
        read_qmelfile(gp,gv,&qcmelt_file,&logPULSEfile);

        // read meteo file
        read_meteofile(gp,gv,&meteo_file,&logPULSEfile);

        // initial conditions
        initiate(gp,gv,&logPULSEfile,&results_flname);

        // call the main PULSE model
        pulsemodel(gp,gv,&logPULSEfile,&results_flname);

        // Simulation completed
        msg = "\n......................................... \n Simulation complete";
        print_screen_log(&logPULSEfile,&msg); 
        
    //} catch(const std::exception& e){
        
    //    try{
    //        msg = "\nError: there was a problem in the code that was not expected (please contact diogo.pinhodacosta@canada.ca)";
    //        print_screen_log(&logPULSEfile,&msg); 
            
    //    }catch(const std::exception& e){
    //    }
      
    //    logPULSEfile.close(); 
        
    //};
    
    return 0;
}

