#include "hydrocalc.h"


/* *****
 * INTERNAL CALCULATION: Determine the addition or removal of upper snow layers  
 * **** */
void watermass_calc_internal(globalvar& gv,globalpar& gp,double* deltt,double *v,
        std::ofstream* logPULSEfile){

    int nl_l = gv.nl;
    int nh_l = gv.nh;
    int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
    //int nt = nli*nhi;
    int ih,il;
    double dv_snow2liqwater,dvfrac_s_dt,existsnow_vol,tofillsnow_vol, 
            add_snow,remove_snow,add_rain,dcomp_swe;

    // timestep fluxes and volumes
    add_snow = std::abs(gv.snowfall_t) * gv.snowl * (*deltt); // as vol mm*mm*m
    remove_snow = std::abs(gv.qmelt_t) * gv.snowl * (*deltt); // as vol mm*mm*m
    add_rain = std::abs(gv.rainfall_t) * gv.snowl * (*deltt); // as vol mm*mm*m

    existsnow_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofillsnow_vol = fmax(gv.v_swe_freshsnow_max - existsnow_vol,0.0f);

    // Refreezing if T<0
    if (gv.tempert_t < 0.0f) 
    {
        (*gv.c_s) = ((*gv.c_s) % (*gv.v_swe) + (*gv.c_m) % (*gv.v_liq))
            / ((*gv.v_swe) + (*gv.v_liq)); // c_m mass will go to c_i
        (*gv.c_m) =  (*gv.c_m) * 0.0f;

        (*gv.v_swe) = (*gv.v_swe) + (*gv.v_liq);
        (*gv.v_liq) = (*gv.v_liq) * 0.0f;

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
    {

        if (tofillsnow_vol >= add_snow) // no need to add_snow new layer
        {
            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                 arma::ones(gv.nl,1) * add_snow * gv.precip_c_t)
                / ((*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add_snow);
            (*gv.v_swe).col(0) = (*gv.v_swe).col(0) + arma::ones(gv.nl,1) * add_snow;
        } 
        else // add_snow new top layer
        {

            (*gv.c_s).col(0) = ((*gv.c_s).col(0) % (*gv.v_swe).col(0) +
                tofillsnow_vol * gv.precip_c_t) / gv.v_swe_freshsnow_max ;
            (*gv.v_swe).col(0) = arma::ones(gv.nl,1) * gv.v_swe_freshsnow_max;
            
            // new layer
            gv.nh++; // remove_snow one layer
            gv.snowH += gv.snowh; // snowpack depth 
            //gv.nh_change -= gv.snowh; 
            gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev + 1,0);
            gv.wetfront_cell = std::max(gv.wetfront_cell + 1,0);

            (*gv.c_m).insert_cols(0,1); // set to zero by default
            (*gv.c_s).insert_cols(0,1); // set to precip_c_t
            (*gv.c_s).col(0) = arma::ones<arma::vec>(nl_l,1) * gv.precip_c_t;

            (*gv.exchange_is).insert_cols(0,1);  // set to zero by default

            (*gv.vfrac2d_m).insert_cols(0,1);  // set to zero by default
            (*gv.vfrac2d_s).insert_cols(0,1); //set to one
            (*gv.vfrac2d_s).col(0) = arma::ones<arma::vec>(nl_l,1);

            (*gv.velc_2d).insert_cols(0,1); 
            (*gv.disp_2d).insert_cols(0,1); 

            (*gv.v_liq).insert_cols(0,1); // set to zero by default
            (*gv.v_swe).insert_cols(0,1);
            (*gv.v_swe).col(0) = arma::ones<arma::vec>(nl_l,1) * (add_snow - tofillsnow_vol);
            (*gv.v_air).insert_cols(0,1);
            (*gv.v_air).col(0) = arma::ones<arma::vec>(nl_l,1) * (gv.vfrac_air_frshsnow * gv.snowl * gv.snowh);
        }

    }

    existsnow_vol = (*gv.v_swe)(round(gv.nl/2)-1,0); 
    tofillsnow_vol = fmax(gv.v_swe_freshsnow_max - existsnow_vol,0.0f);
    
    // Snowmelt or rainfall
    if (remove_snow > 0.0f || add_rain > 0.0f) // melt
    {
        // wetting front calculation
        gv.wetfront_cell_prev = gv.wetfront_cell;
        gv.wetfront_z = std::fmax(gv.wetfront_z - (*v) * (*deltt),0.0f);
        int tmp = nh_l - std::ceil(gv.wetfront_z/gv.snowh);
        gv.wetfront_cell = std::min(tmp,nh_l); // finding the cell when the wetting front is located

        gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev,0);
        gv.wetfront_cell = std::max(gv.wetfront_cell,0);

        if (remove_snow > 0.0f){
            if (existsnow_vol >= remove_snow) // don't remove_snow layer
            {
                

                (*gv.c_m).col(0) = ((*gv.c_m).col(0) % (*gv.v_liq).col(0) + 
                    (*gv.c_s).col(0) * remove_snow)
                    / ( (*gv.v_liq).col(0) + remove_snow);
                // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

                (*gv.v_swe).col(0) = (*gv.v_swe).col(0) - remove_snow;
                (*gv.v_liq).col(0) = (*gv.v_liq).col(0) + remove_snow;
                            
            }else // remove_snow layer
            {
                if ((*gv.c_m).n_cols <= 1) {
                    std::cout << "Snowpack melted completely: sum(qmelt) > sum(precip)" << std::endl;
                }else {
                    // transfer mass from upper layer to lower layer because it is going to be removed
                    (*gv.c_m).col(1) = 
                        ((*gv.c_m).col(0) % (*gv.v_liq).col(0)
                        + (*gv.c_s).col(0) % (*gv.v_swe).col(0) 
                        + (*gv.c_m).col(1) % (*gv.v_liq).col(1)
                        + (*gv.c_s).col(1) % (arma::ones(gv.nl,1) * remove_snow - (*gv.v_swe).col(0))) 
                        / ((*gv.v_liq).col(0) 
                            + (*gv.v_swe).col(0) 
                            + (*gv.v_liq).col(1) 
                            + arma::ones(gv.nl,1) * remove_snow - (*gv.v_swe).col(0));
                    // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

                    (*gv.v_liq).col(1) = (*gv.v_liq).col(1) 
                        + (*gv.v_liq).col(0) 
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
                    (*gv.v_liq).shed_col(0);
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
         (*gv.c_m).col(0) = ((*gv.c_m).col(0) % (*gv.v_liq).col(0) + 
                    add_rain * gv.precip_c_t)
                    / ( (*gv.v_liq).col(0) + add_rain);
                // (*gv.c_s) -> no need to calculate for c_s because it will not change (water masses cancel out)

         (*gv.v_liq).col(0) = (*gv.v_liq).col(0) + add_rain;
        }

        // advection (only water)
        if ((*gv.c_m).n_cols >= 1) {
            for (ih=0;ih<gv.wetfront_cell-1;ih++){
                for (il=0;il<nl_l-1;il++){
                    
                    dv_snow2liqwater = fmin((*v) * (*deltt),1) * (*gv.v_liq).at(il,ih);
                    (*gv.v_liq).at(il,ih) = (*gv.v_liq).at(il,ih) - dv_snow2liqwater;
                    (*gv.v_liq).at(il,ih+1) = (*gv.v_liq).at(il,ih+1) + dv_snow2liqwater;

                }
            }
        }
    }

    // Compaction (Hooke's law)
    if ((*gv.c_m).n_cols > 1){
        for (ih=0;ih<gv.nh-2;ih++){
            for (il=0;il<gv.nl-1;il++){
                
                dcomp_swe = (ih+1) * gv.snowh * gv.compatfact * (gv.v_swe_comp_max - (*gv.v_swe)(il,ih+1));
                dcomp_swe = fmax(0.0f,dcomp_swe);
                dcomp_swe = fmin(dcomp_swe,0.9*(*gv.v_swe)(il,ih));
                
                (*gv.c_s)(il,ih+1) = ((*gv.v_swe)(il,ih+1) * (*gv.c_s)(il,ih+1) 
                                    + dcomp_swe *(*gv.c_s)(il,ih))
                                    / ((*gv.v_swe)(il,ih+1) + dcomp_swe);
                (*gv.c_s)(il,ih) = ((*gv.v_swe)(il,ih) * (*gv.c_s)(il,ih) 
                                    - dcomp_swe * (*gv.c_s)(il,ih))
                                    / ((*gv.v_swe)(il,ih) - dcomp_swe);
                
                (*gv.v_swe)(il,ih) = (*gv.v_swe)(il,ih) - dcomp_swe;
                (*gv.v_swe)(il,ih+1) = (*gv.v_swe)(il,ih+1) + dcomp_swe;
            }
            if ( (*gv.v_swe)(round(gv.nl/2),ih) <= gv.v_swe_comp_min){

                    // transfer mass from upper layer to lower layer because it is going to be removed
                if ((*gv.v_liq)(round(gv.nl/2),ih)  > 0.0f){
                    (*gv.c_m).col(ih+1) = 
                    ((*gv.c_m).col(ih) % (*gv.v_liq).col(ih)
                    + (*gv.c_m).col(ih+1) % (*gv.v_liq).col(ih+1))
                    / ((*gv.v_liq).col(ih) + (*gv.v_liq).col(ih+1));
                }
                
                (*gv.c_s).col(ih+1) = 
                    ((*gv.c_s).col(ih) % (*gv.v_swe).col(ih)
                    + (*gv.c_s).col(ih+1) % (*gv.v_swe).col(ih+1))
                    / ((*gv.v_swe).col(ih) + (*gv.v_swe).col(ih+1));

                (*gv.v_liq).col(ih+1) = (*gv.v_liq).col(ih+1) + (*gv.v_liq).col(ih);
                (*gv.v_swe).col(ih+1) = (*gv.v_swe).col(ih+1) + (*gv.v_swe).col(ih);

                // adust some variables to the removal of a layer
                gv.nh = std::fmax(gv.nh - 1,0);
                gv.wetfront_cell_prev = std::max(gv.wetfront_cell_prev - 1,0);
                gv.wetfront_cell = std::max(gv.wetfront_cell - 1,0);
                gv.snowH = std::fmax(gv.snowH - gv.snowh,0); // snowpack depth

                // remove_snow layer
                (*gv.c_m).shed_col(ih);
                (*gv.c_s).shed_col(ih);
                //(*gv.exchange_si).shed_cols(0,0);
                (*gv.v_liq).shed_col(ih);
                (*gv.v_swe).shed_col(ih);
                (*gv.v_air).shed_col(ih);
                (*gv.exchange_is).shed_col(ih);
                (*gv.vfrac2d_m).shed_col(ih);
                (*gv.vfrac2d_s).shed_col(ih);

                (*gv.velc_2d).shed_col(ih); 
                (*gv.disp_2d).shed_col(ih);

            }
        }
    }
       
    // update mass frac
    (*gv.vfrac2d_m) = (*gv.v_liq) / 
        ((*gv.v_liq) +  (*gv.v_swe) +  (*gv.v_air));
    (*gv.vfrac2d_s) = (*gv.v_swe)
                / ((*gv.v_liq) +  (*gv.v_swe) +  (*gv.v_air));

    // melt volumes and fractions
    if (gv.wetfront_cell > 0)
    {
        gv.vfrac_s = mean(mean((*gv.vfrac2d_s)(arma::span(0,nl_l-1),arma::span(0,gv.wetfront_cell-1))));
        gv.vfrac_m = mean(mean((*gv.vfrac2d_m)(arma::span(0,nl_l-1),arma::span(0,gv.wetfront_cell-1))));
    };


   return;  
}


/* *****
 * EXTERNAL CALCULATION: obtain from external model (e.g., SNOWPACK) -> just supports 1D
 * **** */
bool watermass_calc_external(globalvar& gv,globalpar& gp,double* deltt,
        std::ofstream* logPULSEfile, int t){

    bool err_flag = false;

    // Skip step t = 0 because I need the previous time step
    if (t == 0){
        return err_flag;
    }
    
    //try{
        
        // Get input data at time and t-1
        arma::mat v_swe_ext_t = (*gv.v_swe_ext)(t,arma::span::all);
        arma::mat v_liq_ext_t = (*gv.v_liq_ext)(t,arma::span::all);
        //arma::mat v_ice2liq_1_ext_t = (*gv.v_ice2liq_1_ext)(t,arma::span::all);
        //arma::mat v_ice2liq_2_ext_t = (*gv.v_ice2liq_2_ext)(t,arma::span::all);
        //arma::mat fluxQ_ext_t = (*gv.fluxQ_ext)(t,arma::span::all);
        arma::mat v_swe_ext_tprev = (*gv.v_swe_ext)(t-1,arma::span::all);
        arma::mat v_liq_ext_tprev = (*gv.v_liq_ext)(t-1,arma::span::all);
        //arma::mat v_ice2liq_1_ext_tprev = (*gv.v_ice2liq_1_ext)(t-1,arma::span::all);
        //arma::mat v_ice2liq_2_ext_tprev = (*gv.v_ice2liq_2_ext)(t-1,arma::span::all);
        //arma::mat fluxQ_ext_tprev = (*gv.fluxQ_ext)(t-1,arma::span::all);

        double prec_c_ext_t = (*gv.preci_c_ext)(t,0);

        // Identify the layers with snow
        arma::uvec snowlay_t = arma::find(v_swe_ext_t != 9999);
        arma::uvec snowlay_tprev = arma::find(v_swe_ext_tprev != 9999);
        int num_snowlay_t = snowlay_t.n_elem;
        int num_snowlay_tprev = snowlay_tprev.n_elem;

        // Remove the layers without snow 
        v_swe_ext_t = v_swe_ext_t.cols(snowlay_t);
        v_liq_ext_t = v_liq_ext_t.cols(snowlay_t);
        //v_ice2liq_1_ext_t = v_ice2liq_1_ext_t.cols(snowlay);
        //v_ice2liq_2_ext_t = v_ice2liq_2_ext_t.cols(snowlay);
        //fluxQ_ext_t = fluxQ_ext_t.cols(snowlay);
        v_swe_ext_tprev = v_swe_ext_tprev.cols(snowlay_tprev);
        v_liq_ext_tprev = v_liq_ext_tprev.cols(snowlay_tprev);

        // Determine difference delta_swe and delta_liq (only for the layers in common between t and t-1)
        int num_comon_layers = std::min(num_snowlay_t,num_snowlay_tprev);
        arma::mat delta_swe = arma::reverse(v_swe_ext_t.head_cols(num_comon_layers) - v_swe_ext_tprev.head_cols(num_comon_layers)); // reverse from snowpack to pulse convention
        arma::mat delta_liq = arma::reverse(v_liq_ext_t.head_cols(num_comon_layers) - v_liq_ext_tprev.head_cols(num_comon_layers)); // reverse from snowpack to pulse convention


        // Identify top input (precipitation)
        int diff_num_snowlay = v_swe_ext_t.n_cols - (*gv.v_swe).n_cols;  
        if (diff_num_snowlay > 0){ // Precipitation
            
            gv.nh += diff_num_snowlay; // remove_snow one layer
            gv.snowH += gv.snowh * diff_num_snowlay; // snowpack depth 

            // SWE and liquid
            arma::mat add_snow_layers_swe = arma::reverse(v_swe_ext_t.tail_cols(diff_num_snowlay)); // reverse because pulse adds new layer at col = 0 and not at the end of the array
            arma::mat add_snow_layers_liq = arma::reverse(v_liq_ext_t.tail_cols(diff_num_snowlay)); // reverse because pulse adds new layer at col = 0 and not at the end of the array          
            
            (*gv.vfrac2d_s).insert_cols(0,add_snow_layers_swe); //set to one
            (*gv.vfrac2d_m).insert_cols(0,add_snow_layers_liq);  // set to zero by default
            (*gv.v_swe) = (*gv.vfrac2d_s) * (gv.snowh);
            (*gv.v_liq) = (*gv.vfrac2d_m) * (gv.snowh);
            //(*gv.v_air).insert_cols(0,diff_num_snowlay); // ?? not sure what to put here but likely irrelevatr

            // Other hydraulic variables
            (*gv.exchange_is).insert_cols(0,diff_num_snowlay);  // set to zero by default           

            // Water Quality
            (*gv.c_m).insert_cols(0,diff_num_snowlay); // set to zero by default
            arma::mat add_snow_layers_swe_conc = arma::ones<arma::mat>(1,diff_num_snowlay) * prec_c_ext_t;
            (*gv.c_s).insert_cols(0,add_snow_layers_swe_conc); // set to precip_c_t            
            
        }


        // Mimic mass balance
        double delta_swe_i;
        double delta_liq_i;
        int layer_i;
        for(int il=0;il<num_comon_layers-1 ;il++){

            layer_i = il;
            delta_swe_i = delta_swe(layer_i);
            delta_liq_i = delta_liq(layer_i);

            // delta_swe_i Increase 
            if (delta_swe_i > 0){ // inner loop between top and respective layer to add the necessary masses

                // start by adding to top layer
                //(*gv.c_s).col(0) = ((*gv.c_s).col(0) * (*gv.vfrac2d_s).col(0) + prec_c_ext_t * delta_swe_i) 
                //                / ((*gv.vfrac2d_s).col(0) + delta_swe_i);
                (*gv.vfrac2d_s).col(0) += delta_swe_i;

                // then move that mass from top layer down to the respective layer where an increase was observed
                for (int iil=1;iil<=layer_i ;iil++){ // starting from top now adding the masses

                    //(*gv.c_s).col(iil-1); // concentration here will not change
                    (*gv.vfrac2d_s).col(iil-1) -= delta_swe_i;
                    
                    //(*gv.c_s).col(iil) = ((*gv.c_s).col(iil) * (*gv.vfrac2d_s).col(iil-1) +  delta_swe_i + (*gv.c_s).col(iil-1))
                     //           / ((*gv.vfrac2d_s).col(iil) + delta_swe_i);    
                    (*gv.vfrac2d_s).col(iil) += delta_swe_i;

                }
            }else if(delta_swe_i < 0){  // delta_swe_i Decrease

               
                (*gv.vfrac2d_s).col(layer_i) -= fabs(delta_swe_i);
                (*gv.vfrac2d_m).col(layer_i) += fabs(delta_swe_i);  

            }

            // Liquid water 
            if (delta_liq_i > 0){
                if (layer_i == 0 && delta_swe_i>=0){  // rain

                    (*gv.vfrac2d_m).col(0) += delta_liq_i; // assuming rain because the water needs to come from somewhere

                }else if (layer_i != 0 && delta_swe_i>=0){  // assuming percolation inside the snowpack
                    (*gv.vfrac2d_m).col(layer_i-1) -= delta_liq_i; // remove from upper layer
                    (*gv.vfrac2d_m).col(layer_i) += delta_liq_i; // add to lower layer
                    
                }else if (layer_i != 0 && delta_swe_i<0){ // assuming melt
                    (*gv.vfrac2d_s).col(layer_i) -= delta_liq_i;
                    (*gv.vfrac2d_m).col(layer_i) += delta_liq_i;
                }
            }else if (delta_liq_i < 0 && delta_swe_i >= 0){ // assuming refrezing
                    (*gv.vfrac2d_m).col(layer_i) -= fabs(delta_liq_i);
                    (*gv.vfrac2d_s).col(layer_i) += fabs(delta_liq_i);

            }
        
        }
        (*gv.v_swe) = (*gv.vfrac2d_s) * (gv.snowh);
        (*gv.v_liq) = (*gv.vfrac2d_m) * (gv.snowh);



        // If top layers is removed
        if (diff_num_snowlay < 0){ 
            // Remove the layer
            gv.nh -= abs(diff_num_snowlay);
            gv.snowH -= gv.snowh * abs(diff_num_snowlay);

            (*gv.vfrac2d_m).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.vfrac2d_s).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.v_swe).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.v_liq).shed_cols(0,abs(diff_num_snowlay)-1);
            //(*gv.v_air).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.exchange_is).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.c_m).shed_cols(0,abs(diff_num_snowlay)-1);
            (*gv.c_s).shed_cols(0,abs(diff_num_snowlay)-1);
        }




        /* This calculation was based on v_ice2liq_1_ext_t, v_ice2liq_2_ext_t, adn FluxQ (did not work very well but it might be because of any bugs)
        // First phase change 
        // if ice -> liq (for chem)
        arma::mat v_ice2liq_ext_t_melt = v_ice2liq_1_ext_t; // water mass change
        v_ice2liq_ext_t_melt.elem(arma::find(v_ice2liq_1_ext_t < 0)).zeros(); // only for melt, set freeze to zero here
        arma::mat chemmass_exch_i2l = (*gv.c_s) % v_ice2liq_ext_t_melt * gv.snowh; // calc the chem mass exchange
        arma::mat new_water_mass = ((*gv.vfrac2d_m) + v_ice2liq_ext_t_melt)  * (gv.snowh); // new liquid water mass after applying the phase change
        (*gv.c_m) = ( (*gv.c_m) % (*gv.v_liq) + chemmass_exch_i2l) / new_water_mass; // calculate new chem concentrationj
        arma::uvec non_zero_loc = arma::find(new_water_mass <= gp.num_stblty_thrshld_prsity); // set to zero cases with very low new_wter_mass 
        (*gv.c_m).elem(non_zero_loc).zeros(); // set to zero if water mass is too small

        // if liq -> ice (for chem)
        arma::mat v_ice2liq_ext_t_freeze = v_ice2liq_1_ext_t;
        v_ice2liq_ext_t_freeze.elem(arma::find(v_ice2liq_1_ext_t > 0)).zeros();
        v_ice2liq_ext_t_freeze = - v_ice2liq_ext_t_freeze; 
        arma::mat chemmass_exch_l2i = (*gv.c_m) % v_ice2liq_ext_t_freeze * gv.snowh;    
        new_water_mass = ((*gv.vfrac2d_s) + v_ice2liq_ext_t_freeze)  * (gv.snowh);
        (*gv.c_s) = ( (*gv.c_s) % (*gv.v_swe) + chemmass_exch_l2i) / new_water_mass;
        non_zero_loc = arma::find(new_water_mass <= gp.num_stblty_thrshld_prsity);
        (*gv.c_s).elem(non_zero_loc).zeros(); // set to zero if water mass is too small
        
        // general for hydro
        (*gv.vfrac2d_s) = (*gv.vfrac2d_s) - v_ice2liq_1_ext_t;
        (*gv.vfrac2d_m) = (*gv.vfrac2d_m) + v_ice2liq_1_ext_t;
        (*gv.v_swe) = (*gv.vfrac2d_s) * (gv.snowh);
        (*gv.v_liq) = (*gv.vfrac2d_m) * (gv.snowh);

        // Percolation (just liquid phase)
        // chem
        int num_ele_minus_one = (*gv.vfrac2d_m).n_elem-1;
        if (num_ele_minus_one > 0){
            arma::mat chemmass_exch = (*gv.c_m) % fluxQ_ext_t * gv.snowh;
            (*gv.c_m) = ( (*gv.c_m) % (*gv.v_liq) - chemmass_exch) / ( ((*gv.vfrac2d_m) - fluxQ_ext_t)  * (gv.snowh)); // remove from upper layer
            new_water_mass =  ((*gv.vfrac2d_m).tail_cols(num_ele_minus_one) + 
                            fluxQ_ext_t.head_cols(num_ele_minus_one))  * (gv.snowh);
            (*gv.c_m).tail_cols(num_ele_minus_one) = ( (*gv.c_m).tail_cols(num_ele_minus_one) % (*gv.v_liq).tail_cols(num_ele_minus_one) 
                            + chemmass_exch.head_cols(num_ele_minus_one)) / new_water_mass; // add to lower layer
            non_zero_loc = arma::find(new_water_mass <= gp.num_stblty_thrshld_prsity);
            (*gv.c_m).elem(non_zero_loc).zeros(); // set to zero if water mass is too small
            // hydro
            (*gv.vfrac2d_m) = (*gv.vfrac2d_m) - fluxQ_ext_t; // remove from upper layer
            (*gv.vfrac2d_m).tail_cols(num_ele_minus_one) = (*gv.vfrac2d_m).tail_cols(num_ele_minus_one) + fluxQ_ext_t.head_cols(num_ele_minus_one); // add to lower layer
            (*gv.v_liq) = (*gv.vfrac2d_m) * (gv.snowh);
        };

        // Second phase change
        // if ice -> liq (for chem)
        v_ice2liq_ext_t_melt = v_ice2liq_2_ext_t; // water mass change
        v_ice2liq_ext_t_melt.elem(arma::find(v_ice2liq_2_ext_t > 0)).zeros(); // only for melt, set freeze to zero here
        chemmass_exch_i2l = (*gv.c_s) % v_ice2liq_ext_t_melt * gv.snowh; // calc the chem mass exchange
        new_water_mass = ((*gv.vfrac2d_m) + v_ice2liq_ext_t_melt)  * (gv.snowh); // new liquid water mass after applying the phase change
        (*gv.c_m) = ( (*gv.c_m) % (*gv.v_liq) + chemmass_exch_i2l) / new_water_mass; // calculate new chem concentrationj
        non_zero_loc = arma::find(new_water_mass <= gp.num_stblty_thrshld_prsity); // set to zero cases with very low new_wter_mass 
         (*gv.c_m).elem(non_zero_loc).zeros(); // set to zero if water mass is too small

        // if liq -> ice (fr chem)
        v_ice2liq_ext_t_freeze = v_ice2liq_2_ext_t;
        v_ice2liq_ext_t_freeze.elem(arma::find(v_ice2liq_2_ext_t < 0)).zeros();
        v_ice2liq_ext_t_freeze = -v_ice2liq_ext_t_freeze;
        chemmass_exch_l2i = (*gv.c_m) % v_ice2liq_ext_t_freeze * gv.snowh;    
        new_water_mass = ((*gv.vfrac2d_s) + v_ice2liq_ext_t_freeze)  * (gv.snowh);
        (*gv.c_s) = ( (*gv.c_s) % (*gv.v_swe) + chemmass_exch_l2i) / new_water_mass;
        non_zero_loc = arma::find(new_water_mass <= gp.num_stblty_thrshld_prsity);
        (*gv.c_s).elem(non_zero_loc).zeros(); // set to zero if water mass is too small

        // general for hydro
        (*gv.vfrac2d_s) = (*gv.vfrac2d_s) - v_ice2liq_2_ext_t;
        (*gv.vfrac2d_m) = (*gv.vfrac2d_m) + v_ice2liq_2_ext_t;
        (*gv.v_swe) = (*gv.vfrac2d_s) * (gv.snowh);
        (*gv.v_liq) = (*gv.vfrac2d_m) * (gv.snowh);

*/
/*

    arma::uvec meltloc = find((*gv.v_ice2liq_1_ext) > 0); // cells that melt
    arma::uvec freezeloc = find((*gv.v_ice2liq_1_ext) < 0); // cells that freeze

    // 
    // 1st) CONCENTRATIONS 
    //
    // a) Cells that melt -> meltloc  
    // (*gv.c_s) -> c_s concentration doens't change because it is melting
     (*gv.c_m)(meltloc) = ((*gv.c_m)(meltloc) % (*gv.v_liq)(meltloc) + (*gv.c_s)(meltloc) % (*gv.v_ice2liq_1_ext)(meltloc) * (*deltt) )
        / ((*gv.v_liq)(meltloc) + (*gv.v_ice2liq_1_ext)(meltloc)); 

    // b) Cells that freeze -> freezeloc
    // (*gv.c_m) -> c_m concentration doens't change because it is freezing
    (*gv.c_s)(freezeloc) = ((*gv.c_s)(freezeloc) % (*gv.v_swe)(freezeloc) - (*gv.c_m)(freezeloc) % (*gv.v_ice2liq_1_ext)(freezeloc) * (*deltt) )
    / ((*gv.v_swe)(freezeloc) + (*gv.v_ice2liq_1_ext)(freezeloc));

    //
    // 2nd) update new swe and liqwater masses
    //
    (*gv.v_swe) = (*gv.v_swe) - (*gv.v_ice2liq_1_ext) * (*deltt) ;
    (*gv.v_liq) = (*gv.v_liq) + (*gv.v_ice2liq_1_ext) * (*deltt) ;

    */

    //} catch(const std::exception& e){
    //    err_flag = true;
    //}

   return err_flag;

}