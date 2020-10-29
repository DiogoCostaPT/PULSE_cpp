

#include "readin.h"
#include "toolbox.h"

/* *****
 * Read the simset.pulse file (model set up) 
 * **** */
bool read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,double *h_layer,double *l_layer,
                std::string* qmelt_file,std::string* meteo_file, // if SNOWMODEL = internal 
                std::string* time_file, std::string* v_ice_file,std::string* v_liquid_file, std::string* v_ice2liq_1_file, std::string* v_ice2liq_2_file, std::string* fluxQ_file, // if SNOWMODEL = external
                std::ofstream* logPULSEfile,
                int* n_qmelt_file, int* n_meteo_file, 
                int* n_timExt, int* n_maxLayerExt,
                double *vfrac_air_frshsnow, double *compatfact
               )
{
    
    std::string str, msg, str_hydro_solver, str_snowmodel;
    std::ifstream file(modset_flname);

    bool err_flag = false;
    
    while (std::getline(file, str)) 
    {
        // Continue if comment
        if(str.find("#") != std::string::npos)
            continue;

        // Search for keywords
        if(str.find("COMMNET") != std::string::npos){*sim_purp = str.substr(strlen("COMMNET")+1);}; // comment
        if(str.find("START_TIME") != std::string::npos){(*gp.start_time) = str.substr(strlen("START_TIME")+1);}; // comment
        if(str.find("END_TIME") != std::string::npos){(*gp.end_time) = str.substr(strlen("END_TIME")+1);}; // comment
        if(str.find("PRINT_STEP") != std::string::npos){(gp.print_step) = std::stoi(str.substr(strlen("PRINT_STEP")+1));}; // print time step

        if(str.find("H_LAY_mm") != std::string::npos){(*h_layer) = std::stof(str.substr(strlen("H_LAY_mm")+1))/1000;};  // average roughness height (m)
        if(str.find("L_LAY_mm") != std::string::npos){(*l_layer) = std::stof(str.substr(strlen("L_LAY_mm")+1))/1000;};  // average roughness height (m)

        if(str.find("DENSITY_ICE") != std::string::npos){(gp.rho_ice) = std::stof(str.substr(strlen("DENSITY_ICE")+1));}; // density of ice 
        if(str.find("DENSITY_WATER") != std::string::npos){(gp.rho_water) = std::stof(str.substr(strlen("DENSITY_WATER")+1));}; // density of water
        if(str.find("DENSITY_FRESHSNOW") != std::string::npos){(gp.rho_freshsnow) = std::stof(str.substr(strlen("DENSITY_FRESHSNOW")+1));}; // density of freshsnow
        if(str.find("VFRAC_AIR_FRESHSNOW") != std::string::npos){(*vfrac_air_frshsnow) = std::stof(str.substr(strlen("VFRAC_AIR_FRESHSNOW")+1));};  // volume fraction of air in snow (%)

        if(str.find("A_D") != std::string::npos){(gp.aD) = std::stof(str.substr(strlen("A_D")+1));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("ALPHA_IE") != std::string::npos){(gp.alphaIE) = std::stof(str.substr(strlen("ALPHA_IE")+1));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("COMPFACTOR") != std::string::npos){(*compatfact) = std::stof(str.substr(strlen("COMPFACTOR")+1));}; // compaction factor
        
        // SNOWMODEL
        if(str.find("SNOWMODEL") != std::string::npos){str_snowmodel = removeSpaces(str.substr(strlen("SNOWMODEL")+1));};  // snow model: internal or external
        // if SNOWMODEL = internal
        if(str.find("QMELT_FILE") != std::string::npos){*qmelt_file = removeSpaces(str.substr(strlen("QMELT_FILE")+1));}; // snowmelt file
        if(str.find("METEO_FILE") != std::string::npos){*meteo_file = removeSpaces(str.substr(strlen("METEO_FILE")+1));}; // precipitation and precipitation chemistry file
        // if SNOWMODEL = external
        if(str.find("TIME_FILE") != std::string::npos){*time_file = removeSpaces(str.substr(strlen("TIME_FILE")+1));}; 
        if(str.find("V_ICE_FILE") != std::string::npos){*v_ice_file = removeSpaces(str.substr(strlen("V_ICE_FILE")+1));}; 
        if(str.find("V_LIQUID_FILE") != std::string::npos){*v_liquid_file = removeSpaces(str.substr(strlen("V_LIQUID_FILE")+1));}; 
        if(str.find("V_ICE2LIQ_1_FILE") != std::string::npos){*v_ice2liq_1_file = removeSpaces(str.substr(strlen("V_ICE2LIQ_1_FILE")+1));}; 
        if(str.find("V_ICE2LIQ_2_FILE") != std::string::npos){*v_ice2liq_2_file = removeSpaces(str.substr(strlen("V_ICE2LIQ_2_FILE")+1));};
        if(str.find("FLUXQ_FILE") != std::string::npos){*fluxQ_file = removeSpaces(str.substr(strlen("FLUXQ_FILE")+1));};

        if(str.find("HYDRO_SOLVER") != std::string::npos){str_hydro_solver = removeSpaces(str.substr(strlen("HYDRO_SOLVER")+1));}; // flow solver

    }
    file.close();

    // Checking the model settings
    msg = "> Checking model settings...";
    print_screen_log(logPULSEfile,&msg); 

    msg = "    > SNOWMELT model: " + str_snowmodel;
    print_screen_log(logPULSEfile,&msg); 
    msg = "    > Checking required files for SNOWMELT model = " + str_snowmodel;
    print_screen_log(logPULSEfile,&msg); 

    // Check which SNOWMODEL model and see if corresponding information/files have been provided and can be loaded
    // SNOWMODEL = internal
    if(str_snowmodel.find("internal") != std::string::npos){ 

        gp.snowmodel = 0;

        // check QMELT_FILE
        if(!(*qmelt_file).empty()){
            msg = "        > QMELT_FILE found: " + (*qmelt_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading QMELT_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*qmelt_file),arma::csv_ascii);
            if(flstatus==true){
                *n_qmelt_file = filedata.n_rows;
                msg = "        > QMELT_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > QMELT_FILE: loading failed";
                print_screen_log(logPULSEfile,&msg);
                err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > QMELT_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

         // check METEO_FILE
        if(!(*meteo_file).empty()){
            msg = "        > METEO_FILE found: " + (*meteo_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading METEO_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*meteo_file),arma::csv_ascii);
            if(flstatus==true){
                *n_meteo_file = filedata.n_rows;
                msg = "        > METEO_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > METEO_FILE: loading failed";
                print_screen_log(logPULSEfile,&msg);
                err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > METEO_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

     // SNOWMODEL = external (e.g., SNOWPACK)
    }else if(str_snowmodel.find("external") != std::string::npos){ 

        gp.snowmodel = 1;

        // check TIME_EXT
        if(!(*time_file).empty()){ 
            msg = "        > TIME_FILE found: " + (*time_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading TIME_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*time_file),arma::csv_ascii);
            if(flstatus==true){
                *n_timExt = filedata.n_rows;
                msg = "        > TIME_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > TIME_FILE: loading failed";
                print_screen_log(logPULSEfile,&msg);
                err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > TIME_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }
    

        // check V_ICE_FILE
        if(!(*v_ice_file).empty()){ 
            msg = "        > V_ICE_FILE found: " + (*v_ice_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading V_ICE_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*v_ice_file),arma::csv_ascii);
            if(flstatus==true){
                *n_maxLayerExt = filedata.n_cols;
                msg = "        > V_ICE_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > V_ICE_FILE: loading failed";
                print_screen_log(logPULSEfile,&msg);
                err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > V_ICE_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

         // check V_LIQUID_FILE
        if(!(*v_liquid_file).empty()){
            msg = "        > V_LIQUID_FILE found: " + (*v_liquid_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading V_LIQUID_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*v_liquid_file),arma::csv_ascii);
            if(flstatus==true){
                msg = "        > V_LIQUID_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > V_LIQUID_FILE: loading failed";
                 print_screen_log(logPULSEfile,&msg);
                 err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > V_LIQUID_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

         // check V_ICE2LIQ_1_FILE
        if(!(*v_ice2liq_1_file).empty()){
            msg = "        > V_ICE2LIQ_1_FILE found: " + (*v_ice2liq_1_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading V_ICE2LIQ_1_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*v_ice2liq_1_file),arma::csv_ascii);
            if(flstatus==true){
                msg = "        > V_ICE2LIQ_1_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > V_ICE2LIQ_1_FILE: loading failed";
                 print_screen_log(logPULSEfile,&msg);
                 err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > V_ICE2LIQ_1_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

         // check V_ICE2LIQ_2_FILE
        if(!(*v_ice2liq_2_file).empty()){
            msg = "        > V_ICE2LIQ_2_FILE found: " + (*v_ice2liq_2_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading V_ICE2LIQ_2_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*v_ice2liq_2_file),arma::csv_ascii);
            if(flstatus==true){
                msg = "        > V_ICE2LIQ_2_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > V_ICE2LIQ_2_FILE: loading failed";
                 print_screen_log(logPULSEfile,&msg);
                 err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > V_ICE2LIQ_2_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }

         // check FLUXQ_FILE
        if(!(*fluxQ_file).empty()){
            msg = "        > FLUXQ_FILE found: " + (*fluxQ_file);
            print_screen_log(logPULSEfile,&msg); 

            // Reading FLUXQ_FILE to get size and allocate correct memory in global
            arma::mat filedata; 
            bool flstatus =  filedata.load((*fluxQ_file),arma::csv_ascii);
            if(flstatus==true){
                msg = "        > FLUXQ_FILE: loading successful";
                print_screen_log(logPULSEfile,&msg);
            }else{
                msg = "        > FLUXQ_FILE: loading failed";
                 print_screen_log(logPULSEfile,&msg);
                 err_flag = true;
                return err_flag;
            }

        }else{
            msg = "        > FLUXQ_FILE: not found";
            print_screen_log(logPULSEfile,&msg); 
            err_flag = true;
            return err_flag;
        }
    };  
    
    
    // HYDRO_SOLVER
    msg = "        > HYDRO_SOLVER: " + str_hydro_solver; 
    print_screen_log(logPULSEfile,&msg); 

    if(str_hydro_solver.find("simplified") != std::string::npos){
        gp.hydro_solver = 0;
    }else if(str_hydro_solver.find("CN") != std::string::npos){   
        gp.hydro_solver = 0;
    }else{
        msg = "                > HYDRO_SOLVER: unkown";
        print_screen_log(logPULSEfile,&msg); 
        err_flag = true;
        return err_flag;
    }


    return err_flag;
}    

/* *****
 * Read meteo file
 * **** */
bool read_meteofile(globalpar& gp,globalvar& gv,std::string* meteo_file,
                    std::ofstream* logPULSEfile)
{
    unsigned int a; 
    double tprec=0.0f,snowfall_calc_i=0.0f,temp_calc_i=0.0f,precs_i=0.0f,
            rainfall_calc_i = 0.0f,precip_conc_i = 0.0f;
    std::string msg;
    bool err_flag = false;
     
    arma::mat filedataM; 
    bool flstatusM =  filedataM.load((*meteo_file),arma::csv_ascii);
    int num_cols = filedataM.n_rows-1;

    if(flstatusM == true) {
        for(a=0;a<num_cols;a++){
            tprec = filedataM(a+1,0);  // t prec seconds
            temp_calc_i = filedataM(a+1,1);  // degrees celsius
            rainfall_calc_i = filedataM(a+1,2);  // mm/deltatime
            snowfall_calc_i = filedataM(a+1,3);  // mm/deltatime
            precip_conc_i = filedataM(a+1,4);  // conc of precip
            
            (*gv.meteoall_int).at(a,0) = fabs(tprec);  
            (*gv.meteoall_int).at(a,1) = temp_calc_i;  
            (*gv.meteoall_int).at(a,2) = fabs(rainfall_calc_i)/1000;  
            (*gv.meteoall_int).at(a,3) = fabs(snowfall_calc_i)/1000;  // mm/deltatime -> m/deltatime
            (*gv.meteoall_int).at(a,4) = fabs(precip_conc_i); 
        }
       (gp.Tmeteofile) = tprec; 
       //msg = "Successful loading the file: " + (*meteo_file);
       //print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*meteo_file);   
        print_screen_log(logPULSEfile,&msg);
        err_flag = true;
    } 

    return err_flag;
     
}

/* *****
 * Read snowmelt file
 * **** */
bool read_qmelfile(globalpar& gp,globalvar& gv,std::string* qmelt_file,
    std::ofstream* logPULSEfile)
{
    unsigned int a; 
    double tmelts=0.0f,qcmelt_i;
    std::string msg;

    bool err_flag = false;
      
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load((*qmelt_file),arma::csv_ascii);
    int num_cols = filedataQ.n_rows-1;

    if(flstatusQ == true) {
        for(a=0;a<num_cols;a++){
            tmelts = filedataQ(a+1,0);  // t melt seconds
            qcmelt_i = filedataQ(a+1,1);  // mm/deltatime
            //cmelt_i = filedataQ(a,2);  // value of melt
            (*gv.qcmel_int).at(a,0) = std::fabs(tmelts);  
            (*gv.qcmel_int).at(a,1) = std::fabs(qcmelt_i)/1000; // mm/deltatime -> m/deltatime
            //(*gv.qcmel_int).at(a,2) = cmelt_i; // hh-> sec
             (gp.Tqmeltfile) = tmelts;
        }
       msg = "Successful loading the file: " + (*qmelt_file);
       print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*qmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        err_flag = true;
    } 
        
    return err_flag;
}

