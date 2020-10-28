
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
    
    double H_local,L_local,h_layer,l_layer,vfrac_air_frshsnow,compatfact;
    int nl,nh;
    int n_qmelt_file = 0, n_meteo_file = 0; // SNOWMODEL == internal
    int n_v_ice_file = 0,  n_v_liquid_file = 0, 
        n_v_ice2liq_1_file = 0,  n_v_ice2liq_2_file = 0,  n_fluxQ_file = 0; // SNOWMODEL = external

    bool err_flag = false;

    std::string sim_purp;
    std::string qmelt_file,meteo_file,msg;
    std::string v_ice_file,v_liquid_file,v_ice2liq_1_file,v_ice2liq_2_file,fluxQ_file;
    std::ofstream logPULSEfile ("log.pulse");
    
    std::string modset_flname (argv[1]);
    std::string results_flname (argv[2]);
    std::string signmg;

    // Assign global parameters
    globalpar gp; 
    
    //try{
    
        msg = " ....................................................... \n" 
              " ....................................................... \n" 
              " #####  #   #  #      #####  #####                   \n" 
              " #   #  #   #  #      #      #                       \n" 
              " #####  #   #  #      #####  ###                     \n" 
              " #      #   #  #          #  #                       \n" 
              " #      #####  #####  #####  #####                   \n" 
              " multi-phase multi-layer snowpack chemistry model    \n" 
              " ....................................................... \n"
              " version 2.1: contact diogo.pinhodacosta@canada.ca       \n" 
              " ....................................................... \n";
 
        print_screen_log(&logPULSEfile,&msg);   

        // READ SIMULATIONS SETTINGS
        err_flag = read_simset(gp,modset_flname,
            &sim_purp,&h_layer,&l_layer,
            &qmelt_file,&meteo_file, // if SNOWMODEL = internal
            &v_ice_file,&v_liquid_file,&v_ice2liq_1_file,&v_ice2liq_2_file,&fluxQ_file,  // if SNOWMODEL = external
            &logPULSEfile,
            &n_qmelt_file,&n_meteo_file,
            &n_v_ice_file,&n_v_liquid_file,&n_v_ice2liq_1_file,&n_v_ice2liq_2_file,&n_fluxQ_file,
            &vfrac_air_frshsnow,&compatfact);  

        // TERMINATE MODEL if problem with input data
        if (err_flag == true){
            msg = "Model aborted: problem with input files";   
            print_screen_log(&logPULSEfile,&msg);
            std::abort();
        }

        // CREATE MESH
        checkmesh2(&H_local,&L_local,&h_layer,&l_layer,&nh,&nl,&logPULSEfile,&results_flname);

        // Asign global variables (heap)
        globalvar gv(nh,nl,n_qmelt_file,n_meteo_file); 
        gv.snowH = H_local;
        gv.snowL = L_local;
        gv.snowh = h_layer;
        gv.snowl = l_layer;
        gv.vfrac_air_frshsnow = vfrac_air_frshsnow;
        gv.compatfact = compatfact;

        // Read input files 
        if (gp.snowmodel == 0){ // SNOWMODEL = internal
            
            err_flag = read_qmelfile(gp,gv,&qmelt_file,&logPULSEfile); // read snowmelt input
            if (err_flag == true){
                std::abort();
            }

            err_flag = read_meteofile(gp,gv,&meteo_file,&logPULSEfile); // read meteo file
             if (err_flag == true){
                std::abort();
            }

        }else if(gp.snowmodel == 1){ // SNOWMODEL = external

            read_matrixes_ext(gp,gv,
                &v_ice_file,&v_liquid_file,&v_ice2liq_1_file,&v_ice2liq_2_file,&fluxQ_file);

        }

        

        // Check if input data is consistent
        if (n_qmelt_file!=n_meteo_file){
            if (n_qmelt_file>n_meteo_file){
                signmg = ">";
                gp.Tsim = gp.Tmeteofile;
            }else{
                signmg = "<";
                gp.Tsim = gp.Tqmeltfile;
            }
            msg = "Input timeseries need to have the same size: '" + (qmelt_file) + "' "
                + signmg + " '" + (meteo_file) + "' " + "-> END_TIME reset to min value";   
            print_screen_log(&logPULSEfile,&msg); 
        }else{
            gp.Tsim = gp.Tqmeltfile;
        }

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

