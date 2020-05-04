#include "hydrocalc.h"

/* *****
 * Calculate volume fractions for the different water phases 
 * **** */
void vol_fract_calc(globalpar& gp,globalvar& gv,double *v,double *deltt)
{
    int nli = gv.nl;
    int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
    //int nt = nli*nhi;
    int ih,il;

    // Effect of melt at the top layer (change of mass and fraction_of_mass)
    double dv_snow2liqwater = gv.qmelt_i * (*deltt); // in top layer
    double dvfrac_s_dt = gv.qmelt_i * (*deltt) / (gv.wetfront_cell+1); // top layer

    gv.vfrac_m_prev = gv.vfrac_m;
    gv.vfrac_s_prev = gv.vfrac_s;
    
    gv.vfrac_s = std::fmax(gv.vfrac_s - dvfrac_s_dt, 0.0f);
    gv.vfrac_m = std::fmin(gv.vfrac_m + dvfrac_s_dt, 1.0f);  

    //if (gv.vfrac_s != 0){
    //    gv.vfrac_i = std::fmax(gv.vfrac_i + dvfrac_s_dt * (*deltt) - dvfrac_i_dt * (*deltt), 0.f); 
    //}else{
    //    gv.vfrac_i = std::fmax(gv.vfrac_i - dvfrac_i_dt * (*deltt) , 0.f);     
    //};

    //if(gv.wetfront_cell>0){
    //    (*gv. vfrac2d_m)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell-1)) =  arma::ones(gv.nl,gv.wetfront_cell) * gv.vfrac_m;
    //    (*gv.vfrac2d_s)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell-1)) = arma::ones(gv.nl,gv.wetfront_cell) * gv.vfrac_s;
    //}

    // Volumes of water, swe and air
    if (gv.qmelt_i > 0.0f){

        // top layer (melt)
        (*gv.vfrac2d_m)(arma::span(0,gv.nl-1),0) -= arma::ones(gv.nl,1) * dvfrac_s_dt;
        (*gv.vfrac2d_s)(arma::span(0,gv.nl-1),0) += arma::ones(gv.nl,1) * dvfrac_s_dt;
        (*gv.v_liqwater)(arma::span(0,gv.nl-1),0) += arma::ones(gv.nl,1) * dv_snow2liqwater;
        (*gv.v_swe)(arma::span(0,gv.nl-1),0) -= arma::ones(gv.nl,1) * dv_snow2liqwater;

        // advection (water only)
        for (ih=0;ih<gv.wetfront_cell-1;ih++){
            for (il=0;il<nli-1;il++){
                   
                    dv_snow2liqwater = (*v) * (*deltt) * (*gv.v_liqwater).at(il,ih);
                    (*gv.v_liqwater).at(il,ih) -= dv_snow2liqwater;
                    (*gv.v_liqwater).at(il,ih+1) += dv_snow2liqwater;

                    dvfrac_s_dt = (*v) * (*deltt) * (*gv. vfrac2d_m).at(il,ih);
                    (*gv.vfrac2d_m).at(il,ih) -= dvfrac_s_dt;
                    (*gv.vfrac2d_m).at(il,ih+1) += dvfrac_s_dt;
                    //(*gv.vfrac2d_s)(il,ih) = ;

                    

            }
        }
    }   
       
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
    
     gv.nh_change += (std::abs(gv.precip_i)-std::abs(gv.qmelt_i))*(*deltt); // cell increment

    // Refreezing
    if (gv.qmelt_i==0.0f){ 
        
        (*gv.c_s) += (*gv.c_m)*gv.vfrac_m/gv.vfrac_s; // c_m mass will go to c_i
        (*gv.c_m) *= 0;
        gv.wetfront_cell= 0; // assumes refreezing
        gv.wetfront_cell_prev = 0;
        gv.wetfront_z = gv.snowH;
        gv.vfrac_m=0.008;
        //gv.vfrac_i=0.001;
        gv.vfrac_s= 1 - gv.vfrac_m;// - gv.vfrac_i;
        gv.vfrac_m_prev=gv.vfrac_m;
        //gv.vfrac_i_prev=gv.vfrac_i;
        gv.vfrac_s_prev=gv.vfrac_s;
        gv.upperboundary_cell_prev = 0;

        (*gv.v_swe) += (*gv.v_liqwater);
        (*gv.v_liqwater) *= 0;
        
    }
    
    double v_swe_ref = (*gv.v_swe).at(5,0); // TO BE CHANGED: this will not work if the model us run for truly 2D melt with different snowdepths lateraly

    // layer add or remove
    if (gv.nh_change<0 && std::abs(gv.nh_change)>=v_swe_ref && gv.nh>0){ // MELT - remove layer
        
         // add all immobile and solid slow that melted from the last cell) 
        //(*gv.c_m)(arma::span(0,gv.nl-1),1) = ((*gv.c_m)(arma::span(0,gv.nl-1),0) * gv.vfrac_m_prev
        //    + (*gv.c_s)(arma::span(0,gv.nl-1),0) * gv.vfrac_s_prev
        //    + (*gv.c_i)(arma::span(0,gv.nl-1),0) * gv.vfrac_i_prev ) / gv.vfrac_m; // gv.vfrac_m;

        //(*gv.c_m)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_s)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_i)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        
        gv.nh = std::fmax(gv.nh - 1,0);
        gv.wetfront_cell_prev = gv.wetfront_cell;
        gv.wetfront_cell = std::max(gv.wetfront_cell - 1,0);
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
        (*gv.exchange_is).shed_cols(0,0);
        (*gv.vfrac2d_m).shed_cols(0,0);
        (*gv.vfrac2d_s).shed_cols(0,0);
        (*gv.v_liqwater).shed_cols(0,0);
        (*gv.v_swe).shed_cols(0,0);
        (*gv.v_air).shed_cols(0,0);
            
    } else if (gv.nh_change>0 && std::abs(gv.nh_change)>=gv.v_swe_newlayer){ // adding a new layer

        gv.nh++; // remove one layer
        gv.snowH += gv.snowh; // snowpack depth 
        gv.nh_change -= gv.snowh; 
        gv.wetfront_cell_prev = gv.wetfront_cell;
        gv.wetfront_cell = std::max(gv.wetfront_cell + 1,0);

        (*gv.c_m).insert_cols(0,1);
        (*gv.c_i).insert_cols(0,1);
        (*gv.c_s).insert_cols(0,1);
        arma::mat newsnowlayer;
        newsnowlayer.ones(nl_l,1);
        (*gv.c_s).col(0) = newsnowlayer * gv.precipc_i;
        (*gv.exchange_si).insert_cols(0,1);
        (*gv.exchange_is).insert_cols(0,1);
        (*gv.vfrac2d_m).insert_cols(0,1);
        (*gv.vfrac2d_s).insert_cols(0,1);
        (*gv.v_liqwater).insert_cols(0,1);
        (*gv.v_swe).insert_cols(0,1);
        (*gv.v_swe).col(0) = newsnowlayer * gv.v_swe_newlayer;
        (*gv.v_air).insert_cols(0,1);
    }
   return;  
}
