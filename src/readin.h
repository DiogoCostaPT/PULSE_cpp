

#ifndef READINH_INCLUDED
#define READINH_INCLUDED

#include<armadillo>
#include "global.h"
#include "outwrite.h"

int read_simset(globalpar& gp,const std::string& filename, 
                std::string* sim_purp,int *H_local,int *L_local,
                int *h_layer,int *l_layer,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);

void read_qcmelt(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);


#endif