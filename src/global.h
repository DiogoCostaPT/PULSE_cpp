
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
         hydro_solver; // 0) Crank Nicholson, 1) Forward-time, Central-diff space;
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
  globalvar(size_t nh, size_t nl,size_t n_qcmelt,size_t n_meteoall) {
    this->nh = nh;
    this->nl = nl;
    this->n_qcmelt = n_qcmelt;
    this->n_meteoall = n_meteoall;
   
    c_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    //c_i = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    c_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    // = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    exchange_is = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    qcmel_ts = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_qcmelt,2));
    meteoall_ts = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_meteoall,5));
    velc_2d = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    disp_2d = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

    vfrac2d_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    vfrac2d_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

    v_liqwater = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    v_swe = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    v_air = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));

  }
    
    size_t nh,nl,n_qcmelt,n_meteoall;
  
    std::unique_ptr<arma::Mat<double>> c_m,c_s,qcmel_ts,meteoall_ts,//exchange_si,
                  exchange_is,velc_2d,disp_2d,vfrac2d_m,vfrac2d_s,v_liqwater,v_swe,v_air;
    
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
            qmelt_i = 0.0f,
            tempert_i = 0.0f,
            rainfall_i = 0.0f,
            snowfall_i = 0.0f, 
            precip_c_i = 0.0f,
            v_swe_freshsnow_max = 0.0f,
            v_swe_comp_max = 0.0f,
            v_swe_comp_min = 0.0f,
            vfrac_air_frshsnow = 0.0f,
            compatfact = 0.0f;

    int wetfront_cell = 0,
        wetfront_cell_prev = 0;
                 
};

#endif