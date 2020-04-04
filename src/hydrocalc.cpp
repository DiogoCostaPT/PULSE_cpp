#include "hydrocalc.h"

/* *****
 * Calculate volume fractions for the different water phases 
 * **** */
void vol_fract_calc(globalpar& gp,globalvar& gv,double *deltt)
{
    
    //double dvfrac_s_dt = (*q) / gv.vtotal_check;
    unsigned int il,ih;
    double dvfrac_s_dt = gv.q_i / (gv.wetfront_z / 1000 * gp.rho_frshsnow_init);
    double dvfrac_i_dt = dvfrac_s_dt;

    (*gv.vfrac_m_prev) = (*gv.vfrac_m);
    gv.vfrac_i_prev = gv.vfrac_i;
    (*gv.vfrac_s_prev) = (*gv.vfrac_s);
    
    (*gv.vfrac_s)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell)) = 
        (*gv.vfrac_s)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell)) 
        - dvfrac_s_dt * (*deltt);
    (*gv.vfrac_m)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell)) = 
        (*gv.vfrac_m)(arma::span(0,gv.nl-1),arma::span(0,gv.wetfront_cell)) 
        + dvfrac_s_dt * (*deltt);

    for(il=0;il<gv.nl ;il++){
        for(ih=0;ih<=gv.wetfront_cell;ih++){
            (*gv.vfrac_m)(il,ih) = std::fmax((*gv.vfrac_s)(il,ih),0.0f);
            (*gv.vfrac_s)(il,ih) = std::fmax((*gv.vfrac_s)(il,ih),0.0f);
        }
    }

    if (gv.vfrac_s != 0){
        gv.vfrac_i = std::fmax(gv.vfrac_i + dvfrac_s_dt * 
                            (*deltt) - dvfrac_i_dt * (*deltt), 0.0f); 
    }else{
        gv.vfrac_i = std::fmax(gv.vfrac_i - dvfrac_i_dt * (*deltt) , 0.0f);     
    };
    
       
}

/* *****
 * Calculate the wetting front 
 * **** */
void wetfront_calc(globalpar& gp,globalvar& gv,double *velmax_wtfrt, double *deltt)
{
    int nh_l = gv.nh;
    gv.wetfront_cell_prev = gv.wetfront_cell; // wetting fron cell in the previous time step
    gv.wetfront_z = std::fmax(gv.wetfront_z - (*velmax_wtfrt) * (*deltt),0.0f);
    int tmp = std::floor(nh_l-gv.wetfront_z/gv.snowh);
    gv.wetfront_cell = std::min(tmp,nh_l); // finding the cell when the wetting front is located
    
    if (gv.wetfront_cell > gv.wetfront_cell_prev+1) {
        std::cout << std::to_string(gv.wetfront_cell) << std::endl;
    } 

}


/* *****
 * Determine the addition or removal of upper snow layers  
 * **** */
void upbound_calc(globalvar& gv,globalpar& gp,double* deltt,std::ofstream* logPULSEfile){

    int nl_l = gv.nl;
    
    gv.layer_incrmt += gv.q_i*(*deltt)*gp.rho_m/gp.rho_frshsnow_init; // cell increment
    
    if (gv.q_i>0.0f && gv.layer_incrmt>=gv.snowh){ // MELT - remove layer
        
         // add all immobile and solid slow that melted from the last cell) 
        (*gv.c_m)(arma::span(0,gv.nl-1),1) = ((*gv.c_m)(arma::span(0,gv.nl-1),0) % 
            (*gv.vfrac_m_prev)(arma::span(0,gv.nl-1),0) +
            + (*gv.c_s)(arma::span(0,gv.nl-1),0) % (*gv.vfrac_s_prev)(arma::span(0,gv.nl-1),0)
            + (*gv.c_i)(arma::span(0,gv.nl-1),0) * gv.vfrac_i_prev ) 
            / (*gv.vfrac_m)(arma::span(0,gv.nl-1),0); // gv.vfrac_m;

        //(*gv.c_m)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_s)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        //(*gv.c_i)(arma::span(0,gv.nl-1),arma::span(0,gv.upperboundary_cell_prev)) *= 0;
        
        gv.nh--; // remove one layer
        gv.wetfront_cell--;
        gv.wetfront_cell_prev--;
        gv.snowH -= gv.snowh; // snowpack depth
        gv.wetfront_z -= gv.snowh;
        gv.layer_incrmt -= gv.snowh;
        
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
        (*gv.vfrac_m).shed_cols(0,0);
        (*gv.vfrac_m_prev).shed_cols(0,0);
        (*gv.vfrac_s).shed_cols(0,0);
        (*gv.vfrac_s_prev).shed_cols(0,0);
        (*gv.velc).shed_cols(0,0);
        (*gv.disp).shed_cols(0,0);
        
            
    } else if (gv.q_i<0.0f){ // accumulation period
        
        //(*gv.c_i) += (*gv.c_m)*gv.vfrac_m/gv.vfrac_i; // c_m mass will go to c_i
        (*gv.c_m).reset();
        gv.wetfront_cell= 0; // assumes refreezing
        gv.wetfront_cell_prev = 0;
        gv.wetfront_z = gv.snowH;
        //(*gv.vfrac_m) =0.008;
        //gv.vfrac_i=0.001;
        //gv.vfrac_s= 1 - gv.vfrac_m - gv.vfrac_i;
        //gv.vfrac_m_prev=gv.vfrac_m;
        //gv.vfrac_i_prev=gv.vfrac_i;
        //gv.vfrac_s_prev=gv.vfrac_s;
        gv.upperboundary_cell_prev = 0;
        
        if (abs(gv.layer_incrmt)>gv.snowh){ // adding a new layer

            gv.nh++; // remove one layer
            gv.snowH += gv.snowh; // snowpack depth 
            gv.layer_incrmt += gv.snowh; 

            (*gv.c_m).insert_cols(0,0);
            (*gv.c_i).insert_cols(0,0);
            (*gv.c_s).insert_cols(0,0);
            arma::mat newsnowlayer;
            newsnowlayer.ones(nl_l,1);
            (*gv.c_s).col(0) = newsnowlayer * gv.qc_i;
            (*gv.exchange_si).insert_cols(0,0);
            (*gv.exchange_im).insert_cols(0,0);
            (*gv.vfrac_m).insert_cols(0,0);
            (*gv.vfrac_m_prev).insert_cols(0,0);
            (*gv.vfrac_s).insert_cols(0,0);
            (*gv.vfrac_s_prev).insert_cols(0,0);
            (*gv.velc).insert_cols(0,0);
            (*gv.disp).insert_cols(0,0);
        }
    }  
   return;  
}
