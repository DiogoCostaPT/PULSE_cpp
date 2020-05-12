
#include "initiate.h"

/* *****
 * MAIN  
 * **** */
void initiate(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile,
    std::string* results_flname)
{
    
    unsigned int a, ih, il;
    int nh_l = gv.nh;
    
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
    gv.v_swe_max = (gv.snowh * gv.snowl) * gp.rho_frshsnow_init/gp.rho_m;
    
    if(flstatus == true) 
    {
        for(a=0;a<filedata.col(1).n_elem;a++)
        {
            ih = nh_l - filedata(a,0) - 1;  
            il = filedata(a,1);  
            (*gv.c_m).at(il,ih) = filedata(a,2);
            //(*gv.c_i).at(il,ih) = filedata(a,3);
            (*gv.c_s).at(il,ih) = filedata(a,3);
            (*gv. vfrac2d_m).at(il,ih) = filedata(a,4);
            (*gv. vfrac2d_s).at(il,ih) = filedata(a,5);
            (*gv.v_liqwater).at(il,ih) = filedata(a,6)/(1000*1000); // mm*mm*m -> m*m*m
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

        msg = "Initial conditions found: " + init_file;
        print_screen_log(logPULSEfile,&msg);  
        
    }else{
        msg = "Initial conditions NOT FOUND: simulation aborted";  
        print_screen_log(logPULSEfile,&msg);  
        std::abort();
    }  
    
}