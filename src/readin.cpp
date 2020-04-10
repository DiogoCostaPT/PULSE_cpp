

#include "readin.h"
#include "outwrite.h"

/* *****
 * Read the simset.pulse file (model set up) 
 * **** */
void read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,int *H_local,int *L_local,
                int *h_layer,int *l_layer,std::string* qcmelt_file,
                std::string* meteo_file ,std::ofstream* logPULSEfile,
                int* n_qcmelt, int* n_snowfallt)
{
    
    std::string str, msg, str_hydro_solver;
    
    std::ifstream file(modset_flname);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;

        if(str.find("COMMNET") != std::string::npos){*sim_purp = str.substr(8);}; // comment
        if(str.find("H_LAY") != std::string::npos){(*h_layer) = std::stoi(str.substr(6));};  // average roughness height (m)
        if(str.find("L_LAY") != std::string::npos){(*l_layer) = std::stoi(str.substr(6));};  // average roughness height (m)
        if(str.find("QMELT_FILE") != std::string::npos){*qcmelt_file = str.substr(11);}; // snowmelt file
        if(str.find("METEO_FILE") != std::string::npos){*meteo_file = str.substr(11);}; // snowmelt file
        if(str.find("PRINT_STEP") != std::string::npos){(gp.print_step) = std::stoi(str.substr(11));}; // print time step
        if(str.find("A_D") != std::string::npos){(gp.aD) = std::stof(str.substr(4));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("ALPHA_IE") != std::string::npos){(gp.alphaIE) = std::stof(str.substr(9));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("HYDRO_SOLVER") != std::string::npos){str_hydro_solver = str.substr(13);}; // snowmelt file
    }
    file.close();
    
    if(i==9){
        msg = "Successful loading the file: " + modset_flname;
    } else{
        msg = "PROBLEM loading the file: " + modset_flname;
    } 
    print_screen_log(logPULSEfile,&msg); 
    
    gp.aD /= 3600; // ???
    gp.alphaIE /=3600; // ??

    // Identify the type of hydraulic simulation
    if(str_hydro_solver.find("simplified") != std::string::npos)
        gp.hydro_solver = 0;
    else if(str_hydro_solver.find("snowpack_model") != std::string::npos)
        gp.hydro_solver = 1;

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
        *n_snowfallt = filedata.col(1).n_elem;
    }else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
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
    double tprec=0.0f,prec_i=0.0f,precs_i=0.0f;
    std::string msg;
     
    arma::mat filedataM; 
    bool flstatusM =  filedataM.load((*meteo_file),arma::csv_ascii);
    if(flstatusM == true) {
        for(a=0;a<filedataM.col(1).n_elem;a++){
            tprec = filedataM(a,0);  // t melt seconds
            prec_i = filedataM(a,1);  // value of melt
            precs_i = filedataM(a,2);  // value of melt
            (*gv.snowfall_ts).at(a,0) = tprec;  
            (*gv.snowfall_ts).at(a,1) = prec_i;
            (*gv.snowfall_ts).at(a,2) = precs_i; // hh-> sec
        }
       (gp.Tperd) = tprec;
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
    if(flstatusQ == true) {
        for(a=0;a<filedataQ.col(1).n_elem;a++){
            tmelts = filedataQ(a,0);  // t melt seconds
            qcmelt_i = filedataQ(a,1);  // value of melt
            //cmelt_i = filedataQ(a,2);  // value of melt
            (*gv.qcmel_ts).at(a,0) = tmelts;  
            (*gv.qcmel_ts).at(a,1) = qcmelt_i/3600; // hh-> sec
            //(*gv.qcmel_ts).at(a,2) = cmelt_i; // hh-> sec
        }
       msg = "Successful loading the file: " + (*qcmelt_file);
       print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    } 
        
}

