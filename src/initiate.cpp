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

#include "initiate.h"

/* *****
 * MAIN  
 * **** */
bool initiate(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile,
    std::string* results_flname)
{

    unsigned int a, ih, il;
    int nh_l = gv.nh;
    bool err_flag = false;
    
    arma::mat filedata; 
    std::string init_file, msg;

    char results_flname_char[(*results_flname).size()+1];
    strcpy(results_flname_char,(*results_flname).c_str());
    
    gv.timstart = findLastStep(results_flname_char); // list the results files to get the last time step
    
    init_file = *results_flname + '/' + std::to_string(int(gv.timstart)) + ".txt";
    
    bool flstatus = filedata.load(init_file,arma::csv_ascii);


    //gv.upperboundary_z = gv.nh*gv.snowh;
    gv.wetfront_z = gv.nh*gv.snowh;
    gv.wetfront_cell_prev = 0;
    gv.wetfront_cell = 0;
    //gv.upperboundary_cell = 0;

    // max swe of any new layer added on top of the snowpack due to precipitation
    gv.v_swe_freshsnow_max = (gv.snowh * gv.snowl) * gp.rho_freshsnow/gp.rho_water;
    gv.v_swe_comp_max = (gv.snowh * gv.snowl) * gp.rho_ice/gp.rho_water;
    gv.v_swe_comp_min = 0.001 * gv.snowh * gv.snowl;
    
    msg = " > Checking Initial Conditions (IC)...";
    print_screen_log(logPULSEfile,&msg);  

    if(flstatus == true) 
    {

        msg = "     > IC file found: " + init_file;
        print_screen_log(logPULSEfile,&msg);  
        msg = "     > Loading IC...";
        print_screen_log(logPULSEfile,&msg); 

        try{
            for(a=0;a<filedata.n_rows;a++)
            {
                ih = nh_l - filedata(a,0) - 1;  
                il = filedata(a,1);  
                (*gv.c_m).at(il,ih) = filedata(a,2);
                //(*gv.c_i).at(il,ih) = filedata(a,3);
                (*gv.c_s).at(il,ih) = filedata(a,3);
                (*gv. vfrac2d_m).at(il,ih) = filedata(a,4);
                (*gv. vfrac2d_s).at(il,ih) = filedata(a,5);
                (*gv.v_liq).at(il,ih) = filedata(a,6)/(1000*1000); // mm*mm*m -> m*m*m
                (*gv.v_swe).at(il,ih) = filedata(a,7)/(1000*1000); // mm*mm*m -> m*m*m
                (*gv.v_air).at(il,ih) = filedata(a,8)/(1000*1000); // mm*mm*m -> m*m*m
                //(*gv.exchange_si).at(il,ih) = filedata(a,9);
                //(*gv.exchange_is).at(il,ih) = filedata(a,10);

                if((*gv.c_m).at(il,ih)!=0){
                    gv.wetfront_z = std::fmin((gv.nh - (ih+1)) * gv.snowh,gv.wetfront_z);
                    gv.wetfront_cell = std::min(int(std::round(nh_l-gv.wetfront_z/gv.snowh)),nh_l); // finding the cell when the wetting front is located
                    //gv.upperboundary_z = std::fmax((gv.nh - (ih)) * gv.snowh,gv.upperboundary_z);  
                    //gv.upperboundary_cell = std::min(int(std::round(nh_l-gv.upperboundary_z/gv.snowh)),nh_l);
                }
            }
            
            gv.wetfront_cell_prev = gv.wetfront_cell;

            msg = "     > IC loading: successful";
        
        } catch(const std::exception& e){
             msg = "     > IC loading: failed";
             err_flag = true;
        }
         print_screen_log(logPULSEfile,&msg);  
        
    }else{
        msg = "     > IC file NOT FOUND: " + init_file;
        print_screen_log(logPULSEfile,&msg);  
        err_flag = true;
    }  

    return err_flag;
    
}