

#ifndef READINH_INCLUDED
#define READINH_INCLUDED

#include<armadillo>
#include "global.h"
#include "outwrite.h"

void read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,int *H_local,int *L_local,
                int *h_layer,int *l_layer,std::string* qcmelt_file,
                std::string* meteo_file ,std::ofstream* logPULSEfile,
                int* n_qcmelt, int* n_snowfallt);

void read_qmelfile(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);

void read_meteofile(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);               


#endif