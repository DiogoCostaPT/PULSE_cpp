

#ifndef READINH_INCLUDED
#define READINH_INCLUDED

#include<armadillo>
#include "global.h"
#include "outwrite.h"
#include "toolbox.h"

void read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,double *h_layer,double *l_layer,
                std::string* qcmelt_file,std::string* meteo_file, // if SNOWMODEL = internal 
                std::string* v_ice_file,std::string* v_liquid_file, std::string* v_ice2liq_1, std::string* v_ice2liq_2, // if SNOWMODEL = external
                std::ofstream* logPULSEfile,int* n_qcmelt, 
                int* n_meteoall, double *vfrac_air_frshsnow, double *compatfact
                );

void read_qmelfile(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);

void read_meteofile(globalpar& gp,globalvar& gv,std::string* qcmelt_file,
                std::ofstream* logPULSEfile);               


#endif