// Copyright 2021: Diogo Costa

// This program, PULSE_cpp, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef GLOBALH_INCLUDED
#define GLOBALH_INCLUDED

#include<armadillo>
#include<memory> 

/* ***** 
 * Global Parameters 
 ***** */
class globalpar
{
public:
    
    double Courant=0.8,aD,
           rho_ice=0, // kg.m-3 at 0 degrees
           rho_water=0, // kg.m-3 at 0 degrees
           rho_freshsnow = 320,
           wetfront_z,num_stblty_thrshld_prsity = 1E-6,alphaIE,Tsim,Tmeteofile,Tqmeltfile;
        
    int flag_sens,run_id,s,print_step,
         hydro_solver, // 0) Crank Nicholson, 1) Forward-time, Central-diff space;
         snowmodel; // 0) internal, 1) external
    //std::ofstream logPULSEfile;

    std::string start_time[1],end_time[1];
    
};

/* *****
 * Global Variables 
 ***** */
class globalvar
{
public:
  globalvar() {

  }
  globalvar(size_t nh, size_t nl,size_t n_qmelt_file,size_t n_meteo_file, 
          size_t n_timExt, size_t n_maxLayerExt) {

    this->nh = nh;
    this->nl = nl;
    this->n_qmelt_file = n_qmelt_file;
    this->n_meteo_file = n_meteo_file;
    this->n_timExt = n_timExt;
    this->n_maxLayerExt = n_maxLayerExt;
   
    c_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    //c_i = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    c_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    // = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    exchange_is = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    
    velc_2d = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    disp_2d = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

    vfrac2d_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    vfrac2d_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

    v_liq = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    v_swe = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    v_air = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

    // Depends on SNOWMODEL
    // For snowmodel = internal
    meteoall_int = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_meteo_file,5));
    qcmel_int = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_qmelt_file,2));
    // For snowmodel = external
    time_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,1));
    preci_c_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,1));
    v_liq_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,n_maxLayerExt));
    v_swe_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,n_maxLayerExt));
    v_ice2liq_1_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,n_maxLayerExt));
    v_ice2liq_2_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,n_maxLayerExt));
    fluxQ_ext = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_timExt,n_maxLayerExt));

  }
    
    size_t nh,nl,n_qmelt_file,n_meteo_file,n_timExt, n_maxLayerExt;
  
    std::unique_ptr<arma::Mat<double>> c_m,c_s,exchange_is,velc_2d,disp_2d,vfrac2d_m,vfrac2d_s;
    std::unique_ptr<arma::Mat<double>> v_swe,v_air,v_liq;
    
    std::unique_ptr<arma::Mat<double>> qcmel_int,meteoall_int; // SNOWMODEL = internal
    std::unique_ptr<arma::Mat<double>> time_ext, v_liq_ext, v_swe_ext, 
      v_ice2liq_1_ext, v_ice2liq_2_ext, fluxQ_ext, preci_c_ext; // SNOWMODEL = external

    double snowH = 0.0f, // snowpack depth
            snowL = 0.0f, // snowpack horizontal lenght
            snowl = 0.0f, // grid h lenght
            snowh = 0.0f, // grid l lenght
            vfrac_m= 0.0f,
            vfrac_a = 0.0008f,
            //vfrac_i=0.001,
            vfrac_s= 1 - vfrac_m,// - vfrac_i,
            vfrac_m_prev=vfrac_m,
            //vfrac_i_prev=vfrac_i,
            vfrac_s_prev=vfrac_s,
            timstart = 0.0f,
            wetfront_z = 0.0f,
            //nh_change = 0.0f,
            qmelt_t = 0.0f,
            tempert_t = 0.0f,
            rainfall_t = 0.0f,
            snowfall_t = 0.0f, 
            precip_c_t = 0.0f,
            v_swe_freshsnow_max = 0.0f,
            v_swe_comp_max = 0.0f,
            v_swe_comp_min = 0.0f,
            vfrac_air_frshsnow = 0.0f,
            compatfact = 0.0f;

    int wetfront_cell = 0,
        wetfront_cell_prev = 0;
                 
};

#endif