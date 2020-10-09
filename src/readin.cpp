

#include "readin.h"
#include "outwrite.h"

/* *****
 * Read the simset.pulse file (model set up) 
 * **** */
void read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,double *h_layer,double *l_layer,
                std::string* qcmelt_file,std::string* meteo_file ,
                std::ofstream* logPULSEfile,int* n_qcmelt, 
                int* n_meteoall, double *vfrac_air_frshsnow, double *compatfact
               )
{
    
    std::string str, msg, str_hydro_solver, str_snowmodel;
    
    std::ifstream file(modset_flname);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;

        if(str.find("COMMNET") != std::string::npos){*sim_purp = str.substr(8);}; // comment
        if(str.find("START_TIME") != std::string::npos){(*gp.start_time) = str.substr(11);}; // comment
        if(str.find("END_TIME") != std::string::npos){(*gp.end_time) = str.substr(9);}; // comment
        if(str.find("H_LAY_mm") != std::string::npos){(*h_layer) = std::stof(str.substr(9))/1000;};  // average roughness height (m)
        if(str.find("L_LAY_mm") != std::string::npos){(*l_layer) = std::stof(str.substr(9))/1000;};  // average roughness height (m)
        if(str.find("SNOWMODEL") != std::string::npos){str_snowmodel = str.substr(10);};  // snow model: internal or external
        if(str.find("VFRAC_AIR_FRESHSNOW") != std::string::npos){(*vfrac_air_frshsnow) = std::stof(str.substr(21));};  // volume fraction of air in snow (%)
        if(str.find("QMELT_FILE") != std::string::npos){*qcmelt_file = str.substr(11);}; // snowmelt file
        if(str.find("METEO_FILE") != std::string::npos){*meteo_file = str.substr(11);}; // snowmelt file
        if(str.find("PRINT_STEP") != std::string::npos){(gp.print_step) = std::stoi(str.substr(11));}; // print time step
        if(str.find("A_D") != std::string::npos){(gp.aD) = std::stof(str.substr(4));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("ALPHA_IE") != std::string::npos){(gp.alphaIE) = std::stof(str.substr(9));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("HYDRO_SOLVER") != std::string::npos){str_hydro_solver = str.substr(13);}; // snowmelt file
        if(str.find("COMPFACTOR") != std::string::npos){(*compatfact) = std::stof(str.substr(11));}; // compaction factor
        if(str.find("DENSITY_ICE") != std::string::npos){(gp.rho_ice) = std::stof(str.substr(12));}; // density of ice 
        if(str.find("DENSITY_WATER") != std::string::npos){(gp.rho_water) = std::stof(str.substr(14));}; // density of water
        if(str.find("DENSITY_FRESHSNOW") != std::string::npos){(gp.rho_freshsnow) = std::stof(str.substr(18));}; // density of freshsnow

    }
    file.close();
    
    if(i==16){
        msg = "Successful loading the file: " + modset_flname;
    } else{
        msg = "PROBLEM loading the file: " + modset_flname;
    } 
    print_screen_log(logPULSEfile,&msg); 
    
    //gp.aD /= 3600; // ???
    //gp.alphaIE /=3600; // ??

    // Identify the type of hydraulic simulation
    if(str_hydro_solver.find("simplified") != std::string::npos)
        gp.hydro_solver = 0;
    else if(str_hydro_solver.find("CN") != std::string::npos)
        gp.hydro_solver = 1;

    // Identify the snow model: internal or external
    if(str_snowmodel.find("internal") != std::string::npos)
        gp.snowmodel = 0;
    else if(str_snowmodel.find("external") != std::string::npos)
        gp.snowmodel = 1;

    // Get qmelt file size
    arma::mat filedata; 
    bool flstatus =  filedata.load((*qcmelt_file),arma::csv_ascii);
    if(flstatus==true){
        *n_qcmelt = filedata.col(1).n_elem;
    }else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    }

    // Get meteo file size
    flstatus =  filedata.load((*meteo_file),arma::csv_ascii);
    if(flstatus==true){
        *n_meteoall = filedata.col(1).n_elem;
    }else{
        msg = "PROBLEM loading the file: " + (*meteo_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    }

    return;
}    

/* *****
 * Read meteo file
 * **** */
void read_meteofile(globalpar& gp,globalvar& gv,std::string* meteo_file,
                    std::ofstream* logPULSEfile)
{
    unsigned int a; 
    double tprec=0.0f,snowfall_calc_i=0.0f,temp_calc_i=0.0f,precs_i=0.0f,
            rainfall_calc_i = 0.0f,precip_conc_i = 0.0f;
    std::string msg;
     
    arma::mat filedataM; 
    bool flstatusM =  filedataM.load((*meteo_file),arma::csv_ascii);
    int num_cols = filedataM.col(1).n_elem-1;

    if(flstatusM == true) {
        for(a=0;a<num_cols;a++){
            tprec = filedataM(a+1,0);  // t prec seconds
            temp_calc_i = filedataM(a+1,1);  // degrees celsius
            rainfall_calc_i = filedataM(a+1,2);  // mm/deltatime
            snowfall_calc_i = filedataM(a+1,3);  // mm/deltatime
            precip_conc_i = filedataM(a+1,4);  // conc of precip
            
            (*gv.meteoall_ts).at(a,0) = fabs(tprec);  
            (*gv.meteoall_ts).at(a,1) = temp_calc_i;  
            (*gv.meteoall_ts).at(a,2) = fabs(rainfall_calc_i)/1000;  
            (*gv.meteoall_ts).at(a,3) = fabs(snowfall_calc_i)/1000;  // mm/deltatime -> m/deltatime
            (*gv.meteoall_ts).at(a,4) = fabs(precip_conc_i); 
        }
       (gp.Tmeteofile) = tprec; 
       //msg = "Successful loading the file: " + (*meteo_file);
       //print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*meteo_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    } 
     
}

/* *****
 * Read snowmelt file
 * **** */
void read_qmelfile(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
    std::ofstream* logPULSEfile)
{
    unsigned int a; 
    double tmelts=0.0f,qcmelt_i;
    std::string msg;
      
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load((*qcmelt_file),arma::csv_ascii);
    int num_cols = filedataQ.col(1).n_elem-1;

    if(flstatusQ == true) {
        for(a=0;a<num_cols;a++){
            tmelts = filedataQ(a+1,0);  // t melt seconds
            qcmelt_i = filedataQ(a+1,1);  // mm/deltatime
            //cmelt_i = filedataQ(a,2);  // value of melt
            (*gv.qcmel_ts).at(a,0) = std::fabs(tmelts);  
            (*gv.qcmel_ts).at(a,1) = std::fabs(qcmelt_i)/1000; // mm/deltatime -> m/deltatime
            //(*gv.qcmel_ts).at(a,2) = cmelt_i; // hh-> sec
             (gp.Tqmeltfile) = tmelts;
        }
       msg = "Successful loading the file: " + (*qcmelt_file);
       print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    } 
        
}

