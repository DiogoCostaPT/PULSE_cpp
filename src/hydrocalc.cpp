#include "hydrocalc.h"


/* *****
 * Determine the addition or removal of upper snow layers  
 * **** */
void watermass_calc(globalvar& gv,globalpar& gp,double* deltt,double *v,
        std::ofstream* logPULSEfile){

    int nl_l = gv.nl;
    int nh_l = gv.nh;
    int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
    //int nt = nli*nhi;
    int ih,il;
    double dv_snow2liqwater,dvfrac_s_dt,exist_vol,tofill_vol, add,remove;

    // timestep fluxes and volumes
    add = std::abs(gv.precip_i) * gv.snowl * (*deltt); // as vol mm*mm*m
    remove = std::abs(gv.qmelt_i) * gv.snowl * (*deltt); // as vol mm*mm*m

    exist_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofill_vol = fmax(gv.v_swe_max - exist_vol,0.0f);


    // Precipitation
    if (add > 0.0f)
    {

        if (tofill_vol >= add) // no need to add new layer
        {
            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                 arma::ones(gv.nl,1) * add * gv.precipc_i)
                / ((*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add);
            (*gv.v_swe).col(0) = (*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add;
        } 
        else // add new top layer
        {

            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                tofill_vol * gv.precipc_i) / gv.v_swe_max ;
            (*gv.v_swe).col(0) = arma::ones(gv.nl,1) * gv.v_swe_max;
            
            // new layer
            gv.nh++; // remove one layer
            gv.snowH += gv.snowh; // snowpack depth 
            //gv.nh_change -= gv.snowh; 
            gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev + 1,0);
            gv.wetfront_cell = std::max(gv.wetfront_cell + 1,0);

            (*gv.c_m).insert_cols(0,1); // set to zero by default
            (*gv.c_s).insert_cols(0,1); // set to precipc_i
            (*gv.c_s).col(0) = arma::ones<arma::vec>(nl_l,1) * gv.precipc_i;

            (*gv.exchange_is).insert_cols(0,1);  // set to zero by default

            (*gv.vfrac2d_m).insert_cols(0,1);  // set to zero by default
            (*gv.vfrac2d_s).insert_cols(0,1); //set to one
            (*gv.vfrac2d_s).col(0) = arma::ones<arma::vec>(nl_l,1);

            (*gv.v_liqwater).insert_cols(0,1); // set to zero by default
            (*gv.v_swe).insert_cols(0,1);
            (*gv.v_swe).col(0) = arma::ones<arma::vec>(nl_l,1) * (add - tofill_vol);
            (*gv.v_air).insert_cols(0,1);
            (*gv.v_air).col(0) = arma::ones<arma::vec>(nl_l,1) * (gv.vfrac_air_frshsnow * gv.snowl * gv.snowh);
        }

    }

    exist_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); // need to update because a new layer may have been added
    

    // Melt (and wetfront movement) or Refreezing
    if (remove == 0.0f) // refreezing
    {
        (*gv.c_s) = ((*gv.c_s) % (*gv.v_swe) + (*gv.c_m) % (*gv.v_liqwater))
            / ((*gv.v_swe) + (*gv.v_liqwater)); // c_m mass will go to c_i
        (*gv.c_m) =  (*gv.c_m) * 0.0f;

        (*gv.v_swe) = (*gv.v_swe) + (*gv.v_liqwater);
        (*gv.v_liqwater) = (*gv.v_liqwater) * 0.0f;

        gv.wetfront_cell = 0; // assumes refreezing
        gv.wetfront_cell_prev = 0;
        gv.wetfront_z = gv.snowH;
        gv.vfrac_m = 0.0080f;
        //gv.vfrac_i=0.001;
        gv.vfrac_s = 1.0f - gv.vfrac_m;// - gv.vfrac_i;
        gv.vfrac_m_prev=gv.vfrac_m;
        //gv.vfrac_i_prev=gv.vfrac_i;
        gv.vfrac_s_prev=gv.vfrac_s;
                

    } else // melt
    {

        // wetting front calculation
        gv.wetfront_cell_prev = gv.wetfront_cell;
        gv.wetfront_z = std::fmax(gv.wetfront_z - (*v) * (*deltt),0.0f);
        int tmp = nh_l-std::floor(gv.wetfront_z/gv.snowh);
        gv.wetfront_cell = std::min(tmp,nh_l); // finding the cell when the wetting front is located

        gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev,0);
        gv.wetfront_cell = std::max(gv.wetfront_cell,0);

        if (exist_vol >= remove ) // don't remove layer
        {
            

            (*gv.c_m).col(0) = ((*gv.c_m).col(0) % (*gv.v_liqwater).col(0) + 
                (*gv.c_s).col(0) * remove)
                / ( (*gv.v_liqwater).col(0) + remove);
            // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

            (*gv.v_swe).col(0) = (*gv.v_swe).col(0) - remove;
            (*gv.v_liqwater).col(0) = (*gv.v_liqwater).col(0) + remove;
                        
        }else // remove layer
        {
            if ((*gv.c_m).n_cols <= 1) {
                std::cout << "Snowpack melted completely: sum(qmelt) > sum(precip)" << std::endl;
            }else {
                // transfer mass from upper layer to lower layer because it is going to be removed
                (*gv.c_m).col(1) = 
                    ((*gv.c_m).col(0) % (*gv.v_liqwater).col(0)
                    + (*gv.c_s).col(0) % (*gv.v_swe).col(0) 
                    + (*gv.c_m).col(1) % (*gv.v_liqwater).col(1)
                    + (*gv.c_s).col(1) % (arma::ones(gv.nl,1) * remove - (*gv.v_swe).col(0))) 
                    / ((*gv.v_liqwater).col(0) 
                        + (*gv.v_swe).col(0) 
                        + (*gv.v_liqwater).col(1) 
                        + arma::ones(gv.nl,1) * remove - (*gv.v_swe).col(0));
                // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

                (*gv.v_liqwater).col(1) = (*gv.v_liqwater).col(1) 
                    + (*gv.v_liqwater).col(0) 
                    + (*gv.v_swe).col(0) 
                    + arma::ones(gv.nl,1) * remove - (*gv.v_swe).col(0);
            
                // adust some variables to the removal of a layer
                gv.nh = std::fmax(gv.nh - 1,0);
                gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev - 1,0);
                gv.wetfront_cell = std::max(gv.wetfront_cell - 1,0);
                gv.snowH = std::fmax(gv.snowH - gv.snowh,0); // snowpack depth

                // remove layer
                (*gv.c_m).shed_col(0);
                //(*gv.c_i).shed_cols(0,0);
                (*gv.c_s).shed_col(0);
                //(*gv.exchange_si).shed_cols(0,0);
                (*gv.v_liqwater).shed_col(0);
                (*gv.v_swe).shed_col(0);
                (*gv.v_air).shed_col(0);
                (*gv.exchange_is).shed_col(0);
                (*gv.vfrac2d_m).shed_col(0);
                (*gv.vfrac2d_s).shed_col(0);
                
            }    
        }

        // advection (only water)
        if ((*gv.c_m).n_cols >= 1) {
            for (ih=0;ih<gv.wetfront_cell-1;ih++){
                for (il=0;il<nl_l-1;il++){
                    
                    dv_snow2liqwater = fmin((*v) * (*deltt),1) * (*gv.v_liqwater).at(il,ih);
                    (*gv.v_liqwater).at(il,ih) = (*gv.v_liqwater).at(il,ih) - dv_snow2liqwater;
                    (*gv.v_liqwater).at(il,ih+1) = (*gv.v_liqwater).at(il,ih+1) + dv_snow2liqwater;

                }
            }
        }
    }
    
    // update mass frac
    (*gv.vfrac2d_m) = (*gv.v_liqwater) / 
        ((*gv.v_liqwater) +  (*gv.v_swe) +  (*gv.v_air));
    (*gv.vfrac2d_s) = (*gv.v_swe)
                / ((*gv.v_liqwater) +  (*gv.v_swe) +  (*gv.v_air));

    // melt volumes and fractions
    if (gv.wetfront_cell > 0)
    {
        gv.vfrac_s = mean(mean((*gv.vfrac2d_s)(arma::span(0,nl_l-1),arma::span(0,gv.wetfront_cell-1))));
        gv.vfrac_m = mean(mean((*gv.vfrac2d_m)(arma::span(0,nl_l-1),arma::span(0,gv.wetfront_cell-1))));
    };


   return;  
}
