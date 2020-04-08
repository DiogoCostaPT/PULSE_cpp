

#include "outwrite.h"
#include "readin.h"

/* *****
 * Read the simset.pulse file (model set up) 
 * **** */
int read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,int *H_local,int *L_local,
                int *h_layer,int *l_layer,std::string* qcmelt_file,
                std::ofstream* logPULSEfile)
{
    
    std::string str, msg;
    int n_qcmelt;
    
    std::ifstream file(modset_flname);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;

        if(str.find("COMMNET") != std::string::npos){*sim_purp = str.substr(8);}; // comment
        if(str.find("H_LAY") != std::string::npos){(*h_layer) = std::stoi(str.substr(6));};  // average roughness height (m)
        if(str.find("L_LAY") != std::string::npos){(*l_layer) = std::stoi(str.substr(6));};  // average roughness height (m)
        if(str.find("QMELT_FILE") != std::string::npos){*qcmelt_file = str.substr(11);}; // snowmelt file
        if(str.find("PRINT_STEP") != std::string::npos){(gp.print_step) = std::stoi(str.substr(11));}; // print time step
        if(str.find("A_D") != std::string::npos){(gp.aD) = std::stof(str.substr(4));}; // SWE standard deviation (snow depletion curves, Kevin's paper)
        if(str.find("ALPHA_IE") != std::string::npos){(gp.alphaIE) = std::stof(str.substr(8));}; // SWE standard deviation (snow depletion curves, Kevin's paper)

        gp.aD /= 3600; // ???
        gp.alphaIE /=3600; // ??
    }
    file.close();
    
    if(i==7){
        msg = "Successful loading the file: " + modset_flname;
    } else{
        msg = "PROBLEM loading the file: " + modset_flname;
    } 
    print_screen_log(logPULSEfile,&msg); 
    
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load((*qcmelt_file),arma::csv_ascii);
    if(flstatusQ==true){
        n_qcmelt = filedataQ.col(1).n_elem;
        
    }else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    }
   
    return n_qcmelt;
    
}    

/* *****
 * Read snowmelt and snow-concentration input 
 * **** */
void read_qcmelt(globalpar& gp,globalvar& gv,std::string* qcmelt_file,std::ofstream* logPULSEfile)
{
    unsigned int a; 
    double tmelts=0.0f,tmelts_prev=0.0f,qcmelt_i,cmelt_i;
    std::string msg;
    
    gv.vtotal_check = gv.snowH / 1000 * gp.rho_frshsnow_init; // initial volume
    
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load((*qcmelt_file),arma::csv_ascii);
    if(flstatusQ == true) {
        for(a=0;a<filedataQ.col(1).n_elem;a++){
            tmelts = filedataQ(a,0);  // t melt seconds
            qcmelt_i = filedataQ(a,1);  // value of melt
            cmelt_i = filedataQ(a,2);  // value of melt
            (*gv.qcmelt).at(a,0) = tmelts;  
            (*gv.qcmelt).at(a,1) = qcmelt_i; // hh-> sec
            (*gv.qcmelt).at(a,2) = cmelt_i; // hh-> sec
            gv.vtotal_check += (tmelts-tmelts_prev) * (qcmelt_i/3600); 
            tmelts_prev = tmelts;
        }
       (gp.Tperd) = tmelts;
       msg = "Successful loading the file: " + (*qcmelt_file);
       print_screen_log(logPULSEfile,&msg);
    } else{
        msg = "PROBLEM loading the file: " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    } 
    
    
    if (gv.vtotal_check<0.0f){
        msg = "Snow mass does not balance (initial+accumulation < melt): check the 0.txt file and " + (*qcmelt_file);   
        print_screen_log(logPULSEfile,&msg);
        std::abort();
    }
        
    
    return;
 
}