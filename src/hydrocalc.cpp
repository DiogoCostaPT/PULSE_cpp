#include "hydrocalc.h"

/* *****
 * Calculate volume fractions for the different water phases 
 * **** */
void vol_fract_calc(globalpar& gp,globalvar& gv,double *deltt)
{
    
    //double dvfrac_s_dt = (*q) / gv.vtotal_check;
    double dvfrac_s_dt = gv.qmelt_i / (gv.snowH / 1000 * gp.rho_frshsnow_init);
    double dvfrac_i_dt = dvfrac_s_dt;

    gv.vfrac_m_prev = gv.vfrac_m;
    gv.vfrac_i_prev = gv.vfrac_i;
    gv.vfrac_s_prev = gv.vfrac_s;
    
    gv.vfrac_s = std::fmax(gv.vfrac_s - dvfrac_s_dt * (*deltt), 0.f);
    if (gv.vfrac_s != 0){
        gv.vfrac_i = std::fmax(gv.vfrac_i + dvfrac_s_dt * (*deltt) - dvfrac_i_dt * (*deltt), 0.f); 
    }else{
        gv.vfrac_i = std::fmax(gv.vfrac_i - dvfrac_i_dt * (*deltt) , 0.f);     
    };
    gv.vfrac_m = std::fmin(gv.vfrac_m + dvfrac_i_dt * (*deltt) , 1.f);  
       
}

/* *****
 * Calculate the wetting front 
 * **** */
void wetfront_calc(globalpar& gp,globalvar& gv,double *v, double *deltt)
{
    int nh_l = gv.nh;
    gv.wetfront_cell_prev = gv.wetfront_cell; // wetting fron cell in the previous time step
    gv.wetfront_z = std::fmax(gv.wetfront_z - (*v) * (*deltt),0.0f);
    int tmp = std::round(nh_l-gv.wetfront_z/gv.snowh);
    gv.wetfront_cell = std::min(tmp,nh_l); // finding the cell when the wetting front is located
    
}


/* *****
 * Determine the addition or removal of upper snow layers  
 * **** */
void upbound_calc(globalvar& gv,globalpar& gp,double* deltt,std::ofstream* logPULSEfile){

    int nl_l = gv.nl;
    
    gv.nh_change += (std::abs(gv.precip_i)-std::abs(gv.qmelt_i))*(*deltt)
                    *gp.rho_m/gp.rho_frshsnow_init; // cell increment

    // Refreezing
    if (gv.qmelt_i==0.0f){ 
        
        (*gv.c_i) += (*gv.c_m)*gv.vfrac_m/gv.vfrac_i; // c_m mass will go to c_i
        (*gv.c_m) *= 0;
        gv.wetfront_cell= 0; // assumes refreezing
        gv.wetfront_cell_prev = 0;
        gv.wetfront_z = gv.snowH;
        gv.vfrac_m=0.008;
        gv.vfrac_i=0.001;
        gv.vfrac_s= 1 - gv.vfrac_m - gv.vfrac_i;
        gv.vfrac_m_prev=gv.vfrac_m;
        gv.vfrac_i_prev=gv.vfrac_i;
        gv.vfrac_s_prev=gv.vfrac_s;
        gv.upperboundary_cell_prev = 0;
    }
    
    // layer add or remove
    if (gv.nh_change<0 && std::abs(gv.nh_change)>=gv.snowh){ // MELT - remove layer
        
         // add all immobile and solid slow that melted from the last cell) 
        //(*gv.c_m)(arma::span(0,gv.nl-1),1) = ((*gv.c_m)(arma::span(0,gv.nl-1),0) * gv.vfrac_m_prev
        //    + (*gv.c_s)(arma::span(0,gv.nl-1),0) * gv.vfrac_s_prev
        //    + (*gv.c_i)(arma::span(0,gv.nl-1),0) * gv.vfrac_i_prev ) / gv.vfrac_m; // gv.vfrac_m;

        //(*gv.c_m)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_s)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_i)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        
        gv.nh = std::fmax(gv.nh - 1,0);
        gv.wetfront_cell = std::max(gv.wetfront_cell - 1,0);
        gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev - 1,0);
        gv.snowH = std::fmax(gv.snowH - gv.snowh,0); // snowpack depth
        //gv.wetfront_z -= gv.snowh; -> it should not decrease because this is counted from the bottom
        gv.nh_change += gv.snowh;
        
        //gv.upperboundary_cell_prev = gv.upperboundary_cell; // wetting fron cell in the previous time step
        //gv.upperboundary_z =  std::fmax(gv.upperboundary_z - (*q) * (*deltt),0.0f);
        //tmp_int = int(std::round(gv.nh-gv.upperboundary_z/gv.snowh));
        //gv.upperboundary_cell = std::min(tmp_int,nh_l); // in what cell is the wetting front
        //upperboundary_cell = [upperboundary_cell, upperboundary_cell_new];
        
        (*gv.c_m).shed_cols(0,0);
        (*gv.c_i).shed_cols(0,0);
        (*gv.c_s).shed_cols(0,0);
        (*gv.exchange_si).shed_cols(0,0);
        (*gv.exchange_im).shed_cols(0,0);
            
    } else if (gv.nh_change>0 && std::abs(gv.nh_change)>=gv.snowh){ // adding a new layer

            gv.nh++; // remove one layer
            gv.snowH += gv.snowh; // snowpack depth 
            gv.nh_change -= gv.snowh; 

            (*gv.c_m).insert_cols(0,1);
            (*gv.c_i).insert_cols(0,1);
            (*gv.c_s).insert_cols(0,1);
            arma::mat newsnowlayer;
            newsnowlayer.ones(nl_l,1);
            (*gv.c_s).col(0) = newsnowlayer * gv.precipc_i;
            (*gv.exchange_si).insert_cols(0,1);
            (*gv.exchange_im).insert_cols(0,1);
    }
   return;  
}
