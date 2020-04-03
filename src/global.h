
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
           rho_s=917.0, // kg.m-3 at 0 degrees
           rho_m=998.8, // kg.m-3 at 0 degrees
           rho_frshsnow_init = 320,
           wetfront_z,num_stblty_thrshld_prsity = 1E-6,alphaIE,Tperd;
    
    int flag_sens,run_id,s,print_step;
    //std::ofstream logPULSEfile;
    
};

/* *****
 * Global Variables 
 ***** */
class globalvar
{
public:
  globalvar() {

  }
  globalvar(size_t nh, size_t nl,size_t n_qcmelt) {
    this->nh = nh;
    this->nl = nl;
    this->n_qcmelt = n_qcmelt;
   
    c_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    c_i = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    c_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    exchange_si = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    exchange_im = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nl,nh));
    qcmelt = std::unique_ptr<arma::Mat<double>>( new  arma::mat(n_qcmelt,3));

  }
    
    size_t nh,nl,n_qcmelt,
            snowH, // snowpack depth
            snowL, // snowpack horizontal lenght
            snowl, // grid h lenght
            snowh; // grtid l lenght
  
    std::unique_ptr<arma::Mat<double>> c_m,c_i,c_s,qcmelt,exchange_si,exchange_im;
    double vfrac_m=0.008,
            vfrac_i=0.001,
            vfrac_s= 1 - vfrac_m - vfrac_i,
            vfrac_m_prev=vfrac_m,
            vfrac_i_prev=vfrac_i,
            vfrac_s_prev=vfrac_s,
            vtotal_check,
            timstart,
            wetfront_z,
            wetfront_cell,
            wetfront_cell_prev,
            upperboundary_cell_prev,
            layer_incrmt,q_i,qc_i;
                 
};

#endif