

#ifndef READINH_INCLUDED
#define READINH_INCLUDED

#include<armadillo>
#include "global.h"
#include "outwrite.h"
#include "toolbox.h"

bool read_simset(globalpar& gp,const std::string& modset_flname, 
                std::string* sim_purp,double *h_layer,double *l_layer,
                std::string* qmelt_file,std::string* meteo_file, // if SNOWMODEL = internal 
                std::string* time_ext_file, std::string* v_ice_file,std::string* v_liquid_file, std::string* v_ice2liq_1_file, std::string* v_ice2liq_2_file, std::string* fluxQ_file, // if SNOWMODEL = external
                std::ofstream* logPULSEfile,
                int* n_qmelt_file, int* n_meteo_file, 
                int* n_timExt, int* n_maxLayerExt,
                double *vfrac_air_frshsnow, double *compatfact
               );

bool read_qmelfile(globalpar& gp,globalvar& gv,std::string* qmelt_file,
                std::ofstream* logPULSEfile);

bool read_meteofile(globalpar& gp,globalvar& gv,std::string* qmelt_file,
                std::ofstream* logPULSEfile);               


#endif