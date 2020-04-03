

#include "outwrite.h"
#include "readin.h"

/* *****
 * Read the simset.pulse file (model set up) 
 * **** */
int read_simset(globalpar& gp,std::string* sim_purp, int *H_local,int *L_local, int *h_layer,int *l_layer, std::string* qcmelt_file,std::ofstream* logPULSEfile)
{
    
    std::string str, modset_flname, msg;
    int n_qcmelt;
    modset_flname = "bin/simset.pulse";
    
    std::ifstream file(modset_flname);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;
        if(i==1){(*sim_purp) = str;};
        //if(i==2){(*H_local) = std::round(std::stoi(str));};
        //if(i==3){(*L_local) = std::round(std::stoi(str));};
        if(i==2){(*h_layer) = std::round(std::stoi(str));};
        if(i==3){(*l_layer) = std::round(std::stoi(str));};
        if(i==4){(*qcmelt_file) = str;};  
        if(i==5){gp.print_step = std::stoi(str);};
        if(i==6){gp.aD = std::stof(str)/3600;}; 
        if(i==7){gp.alphaIE = std::stof(str)/3600;}; 
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
            (*gv.qcmelt).at(a,1) = qcmelt_i/3600; // hh-> sec
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