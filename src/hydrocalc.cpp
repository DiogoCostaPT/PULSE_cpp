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
<<<<<<< Updated upstream
    double dv_snow2liqwater,dvfrac_s_dt,existsnow_vol,tofillsnow_vol, 
            add_snow,remove_snow,add_rain,dcomp_swe;
=======
    double dv_snow2liqwater,dvfrac_s_dt,exist_vol,tofill_vol, add,remove;

    // timestep fluxes and volumes
    add = std::abs(gv.precip_i) * gv.snowl * (*deltt); // as vol mm*mm*m
    remove = std::abs(gv.qmelt_i) * gv.snowl * (*deltt); // as vol mm*mm*m
>>>>>>> Stashed changes

    // timestep fluxes and volumes
    add_snow = std::abs(gv.snowfall_i) * gv.snowl * (*deltt); // as vol mm*mm*m
    remove_snow = std::abs(gv.qmelt_i) * gv.snowl * (*deltt); // as vol mm*mm*m
    add_rain = std::abs(gv.rainfall_i) * gv.snowl * (*deltt); // as vol mm*mm*m

<<<<<<< Updated upstream
    existsnow_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofillsnow_vol = fmax(gv.v_swe_freshsnow_max - existsnow_vol,0.0f);

    // Refreezing if T<0
    if (gv.tempert_i < 0.0f) 
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
                

    }

    existsnow_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofillsnow_vol = fmax(gv.v_swe_freshsnow_max - existsnow_vol,0.0f);

    // Snowfall
    if (add_snow > 0.0f)
=======

    // Precipitation
    if (add > 0.0f)
>>>>>>> Stashed changes
    {

        if (tofillsnow_vol >= add_snow) // no need to add_snow new layer
        {
            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                 arma::ones(gv.nl,1) * add_snow * gv.precip_c_i)
                / ((*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add_snow);
            (*gv.v_swe).col(0) = (*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add_snow;
        } 
        else // add_snow new top layer
        {

            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                tofillsnow_vol * gv.precip_c_i) / gv.v_swe_freshsnow_max ;
            (*gv.v_swe).col(0) = arma::ones(gv.nl,1) * gv.v_swe_freshsnow_max;
            
            // new layer
            gv.nh++; // remove_snow one layer
            gv.snowH += gv.snowh; // snowpack depth 
            //gv.nh_change -= gv.snowh; 
            gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev + 1,0);
            gv.wetfront_cell = std::max(gv.wetfront_cell + 1,0);

            (*gv.c_m).insert_cols(0,1); // set to zero by default
            (*gv.c_s).insert_cols(0,1); // set to precip_c_i
            (*gv.c_s).col(0) = arma::ones<arma::vec>(nl_l,1) * gv.precip_c_i;

            (*gv.exchange_is).insert_cols(0,1);  // set to zero by default

            (*gv.vfrac2d_m).insert_cols(0,1);  // set to zero by default
            (*gv.vfrac2d_s).insert_cols(0,1); //set to one
            (*gv.vfrac2d_s).col(0) = arma::ones<arma::vec>(nl_l,1);

            (*gv.velc_2d).insert_cols(0,1); 
            (*gv.disp_2d).insert_cols(0,1); 

            (*gv.v_liqwater).insert_cols(0,1); // set to zero by default
            (*gv.v_swe).insert_cols(0,1);
            (*gv.v_swe).col(0) = arma::ones<arma::vec>(nl_l,1) * (add_snow - tofillsnow_vol);
            (*gv.v_air).insert_cols(0,1);
            (*gv.v_air).col(0) = arma::ones<arma::vec>(nl_l,1) * (gv.vfrac_air_frshsnow * gv.snowl * gv.snowh);
        }

    }

    existsnow_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofillsnow_vol = fmax(gv.v_swe_freshsnow_max - existsnow_vol,0.0f);
    
<<<<<<< Updated upstream
    // Snowmelt or rainfall
    if (remove_snow > 0.0f || add_rain > 0.0f) // melt
=======

    // Melt (and wetfront movement) or Refreezing
    if (remove == 0.0f) // refreezing
    {
        (*gv.c_s) = ( (*gv.c_s) % (*gv.v_swe) + (*gv.c_m) % (*gv.v_liqwater))
            / ((*gv.v_swe) + (*gv.v_liqwater)); // c_m mass will go to c_i
        (*gv.c_m) *= 0;

        (*gv.v_swe) += (*gv.v_liqwater);
        (*gv.v_liqwater) *= 0;

        gv.wetfront_cell = 0; // assumes refreezing
        gv.wetfront_cell_prev = 0;
        gv.wetfront_z = gv.snowH;
        gv.vfrac_m = 0.008;
        //gv.vfrac_i=0.001;
        gv.vfrac_s = 1 - gv.vfrac_m;// - gv.vfrac_i;
        gv.vfrac_m_prev=gv.vfrac_m;
        //gv.vfrac_i_prev=gv.vfrac_i;
        gv.vfrac_s_prev=gv.vfrac_s;
                

    } else // melt
>>>>>>> Stashed changes
    {
        // wetting front calculation
        gv.wetfront_cell_prev = gv.wetfront_cell;
        gv.wetfront_z = std::fmax(gv.wetfront_z - (*v) * (*deltt),0.0f);
        int tmp = nh_l - std::ceil(gv.wetfront_z/gv.snowh);
        gv.wetfront_cell = std::min(tmp,nh_l); // finding the cell when the wetting front is located

        gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev,0);
        gv.wetfront_cell = std::max(gv.wetfront_cell,0);

        if (remove_snow > 0.0f){
            if (existsnow_vol >= remove_snow ) // don't remove_snow layer
            {
                

                (*gv.c_m).col(0) = ((*gv.c_m).col(0) % (*gv.v_liqwater).col(0) + 
                    (*gv.c_s).col(0) * remove_snow)
                    / ( (*gv.v_liqwater).col(0) + remove_snow);
                // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

                (*gv.v_swe).col(0) = (*gv.v_swe).col(0) - remove_snow;
                (*gv.v_liqwater).col(0) = (*gv.v_liqwater).col(0) + remove_snow;
                            
            }else // remove_snow layer
            {
                if ((*gv.c_m).n_cols <= 1) {
                    std::cout << "Snowpack melted completely: sum(qmelt) > sum(precip)" << std::endl;
                }else {
                    // transfer mass from upper layer to lower layer because it is going to be removed
                    (*gv.c_m).col(1) = 
                        ((*gv.c_m).col(0) % (*gv.v_liqwater).col(0)
                        + (*gv.c_s).col(0) % (*gv.v_swe).col(0) 
                        + (*gv.c_m).col(1) % (*gv.v_liqwater).col(1)
                        + (*gv.c_s).col(1) % (arma::ones(gv.nl,1) * remove_snow - (*gv.v_swe).col(0))) 
                        / ((*gv.v_liqwater).col(0) 
                            + (*gv.v_swe).col(0) 
                            + (*gv.v_liqwater).col(1) 
                            + arma::ones(gv.nl,1) * remove_snow - (*gv.v_swe).col(0));
                    // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

                    (*gv.v_liqwater).col(1) = (*gv.v_liqwater).col(1) 
                        + (*gv.v_liqwater).col(0) 
                        + (*gv.v_swe).col(0) 
                        + arma::ones(gv.nl,1) * remove_snow - (*gv.v_swe).col(0);
                
                    // adust some variables to the removal of a layer
                    gv.nh = std::fmax(gv.nh - 1,0);
                    gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev - 1,0);
                    gv.wetfront_cell = std::max(gv.wetfront_cell - 1,0);
                    gv.snowH = std::fmax(gv.snowH - gv.snowh,0); // snowpack depth

                    // remove_snow layer
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

                    (*gv.velc_2d).shed_col(0); 
                    (*gv.disp_2d).shed_col(0);
                    
                }    
            }
        }
        if (add_rain > 0.0f){
         (*gv.c_m).col(0) = ((*gv.c_m).col(0) % (*gv.v_liqwater).col(0) + 
                    add_rain * gv.precip_c_i)
                    / ( (*gv.v_liqwater).col(0) + add_rain);
                // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

         (*gv.v_liqwater).col(0) = (*gv.v_liqwater).col(0) + add_rain;
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

    // Compaction (Hooke's law)
    for (ih=0;ih<nh_l-2;ih++){
            for (il=0;il<nl_l-1;il++){
                
                dcomp_swe = (ih+1) * gv.snowh * gv.compatfact * (gv.v_swe_comp_max - (*gv.v_swe)(il,ih+1));
                dcomp_swe = fmax(0.0f,dcomp_swe);
                dcomp_swe = fmin(dcomp_swe,(*gv.v_swe)(il,ih));
                
                (*gv.c_s)(il,ih+1) = ((*gv.v_swe)(il,ih+1) * (*gv.c_s)(il,ih+1) 
                                    + dcomp_swe *(*gv.c_s)(il,ih))
                                    / ((*gv.v_swe)(il,ih+1) + dcomp_swe);
                (*gv.c_s)(il,ih) = ((*gv.v_swe)(il,ih) * (*gv.c_s)(il,ih) 
                                    - dcomp_swe * (*gv.c_s)(il,ih))
                                    / ((*gv.v_swe)(il,ih) - dcomp_swe);
                
                (*gv.v_swe)(il,ih) = (*gv.v_swe)(il,ih) - dcomp_swe;
                (*gv.v_swe)(il,ih+1) = (*gv.v_swe)(il,ih+1) + dcomp_swe;

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
