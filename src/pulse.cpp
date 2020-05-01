
#include "pulse.h"


/* *****
 * Main PULSE model  
 * **** */
void pulsemodel(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile,
        std::string* results_flname)
{
    
    // initiation
    double tcum = 0.f,
            //q = 0.f, // melt volume/int
            t = 1.f, 
            deltt = 1.0f, // time step calculated from the CFL condition
            velc = 0.0f, // interstitial flow velocity [m s-1]
            D = 0.0f; // dispersion coefficient [m2/s]
    std::string msg;  
    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    bool outwritestatus;
    arma::mat exchange_i;
    
    unsigned int il,ih,print_next;
      
    tcum = gv.timstart;
    print_next = tcum + gp.print_step;
        
    while (tcum < gp.Tperd)
    {

        t += 1;
                
        findInterpPrec(gv,&tcum);
        findInterpQmelt(gv,&tcum);

        if (gv.qmelt_i==0.0f){ // accumulation only 
            deltt = std::fmin(gp.print_step,
                    gv.snowh * gp.rho_frshsnow_init / (std::abs(gv.precip_i)*gp.rho_m)); 
            upbound_calc(gv,gp,&deltt,logPULSEfile);
        } else {// melt       
                // Estimate interstitial flow velocity 
            velc = gv.qmelt_i / gv.vfrac_m; // interstitial flow velocity [m s-1]
            D = gp.aD * velc;       // dispersion coefficient [m2/s]

            (*gv.velc_2d) = (*gv.velc_2d)*0 + velc;
            (*gv.disp_2d) = (*gv.disp_2d)*0 + D;
            
            deltt = std::fmin(gp.Courant * gv.snowh / velc,gp.print_step);

            // limit step so that if there is melt or accumulation it doesn't go more than one cell
            deltt = std::fmin(deltt,
                    gv.snowh * gp.rho_frshsnow_init / (std::abs((std::abs(gv.precip_i)-std::abs(gv.qmelt_i)))*gp.rho_m));
            
            upbound_calc(gv,gp,&deltt,logPULSEfile);
            
        }

        tcum = tcum + deltt; 

        if (gv.qmelt_i>0.0f && gv.nh>0){ // if melt
            
            // wetting front
            wetfront_calc(gp,gv,&velc,&deltt);
            
            // calculate volume fractions
            vol_fract_calc(gp,gv,&velc,&deltt);

            // check CFC validation
            if((gv.wetfront_cell - gv.wetfront_cell_prev) > 1){
                msg = "CFC condition violation - check code";
                print_screen_log(logPULSEfile,&msg);
                abort();
            }

            // 
            //if (gv.vfrac_m < 1-gp.num_stblty_thrshld_prsity && gv.vfrac_i > gp.num_stblty_thrshld_prsity && gv.vfrac_s > gp.num_stblty_thrshld_prsity){
            if (gv.vfrac_m < 1-gp.num_stblty_thrshld_prsity && gv.vfrac_s > gp.num_stblty_thrshld_prsity){  if (gv.wetfront_cell > 5){ // to have sufficient layers for ADE solver

                if (gp.hydro_solver == 0){
                   crank_nicholson(gv,&deltt,&velc,&D); // solve advection and dispersion in the mobile zone
                }else if (gp.hydro_solver == 1){
                   //crank_nicholson_hydr2D(gv,&deltt);
                   FtCs_solve_hydr2D(gv,&deltt);

                    // Crank Nicolson to limit the fluxes across boundaries
                   //exchange_i = arma::max(v * deltt * ((*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1) - (*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new))/gv.snowh,(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1));
                   //(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1) -= exchange_i; // compute onh advection to the wetting front
                   //(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new) += exchange_i; // compute onh advection to the wetting front
                }
              }

            if (gv.qmelt_i > 0.0f){

                 // Ion exclusion: Exchange with immobile phase (just exchange)
                (*gv.exchange_is)  = deltt * (gp.alphaIE * (*gv.c_s) * gv.vfrac_s) ; 

                // limit the flux to the available material
                for(il=0;il<gv.nl ;il++){
                    for(ih=0;ih<gv.nh ;ih++){
                         if ((*gv.exchange_is).at(il,ih) > 0 && ih<gv.wetfront_cell){
                            (*gv.exchange_is).at(il,ih) = std::fmin((*gv.exchange_is).at(il,ih),
                                (*gv.c_s).at(il,ih) * gv.vfrac_s);

                            (*gv.c_s).at(il,ih) = ((*gv.c_s).at(il,ih) * gv.vfrac_s 
                                    - (*gv.exchange_is).at(il,ih)) / gv.vfrac_s;
                            (*gv.c_m).at(il,ih) = ((*gv.c_m).at(il,ih) * gv.vfrac_m 
                                    + (*gv.exchange_is).at(il,ih)) / gv.vfrac_m;
                            (*gv.c_s).at(il,ih) = std::fmax((*gv.c_s).at(il,ih),0.0f);
                            (*gv.c_m).at(il,ih) = std::fmax((*gv.c_m).at(il,ih),0.0f);

                        //}else if((*gv.exchange_is).at(il,ih) < 0){
                        //     (*gv.exchange_is).at(il,ih) = 0.0f;
                             //msg = "PROBLEM: negative s->i exchange" + std::to_string((*gv.c_s).at(il,ih));
                             //std::cout << msg << std::endl;
                             //print_screen_log(logPULSEfile,&msg);
                        }
                         //if(ih<gv.upperboundary_cell || ih>gv.wetfront_cell){
                        //if(ih>gv.wetfront_cell){
                        //    (*gv.exchange_is).at(il,ih) = 0.0f;
                        //};
                    }
                };
                
                //(*gv.c_s) =  ( (*gv.c_s) * gv.vfrac_s - (*gv.exchange_is)) / gv.vfrac_s; // / vfrac_m(t);
                //(*gv.c_m) = ( (*gv.c_m) * gv.vfrac_m + (*gv.exchange_is)) / gv.vfrac_m;
  

                
                // Snowmelt
                (*gv.exchange_si) = (*gv.c_s) * gv.qmelt_i * deltt / (gv.wetfront_cell+1);

                // limit the flux to the available material
                for(il=0;il<gv.nl ;il++){
                    for(ih=0;ih<gv.nh ;ih++){

                            if ((*gv.exchange_si).at(il,ih) > 0 && ih<gv.wetfront_cell){
                                (*gv.exchange_si).at(il,ih) = std::fmin((*gv.exchange_si).at(il,ih),(*gv.c_s).at(il,ih));
                            
                                (*gv.c_s).at(il,ih) = ((*gv.c_s).at(il,ih) * gv.vfrac_s 
                                    - (*gv.exchange_si).at(il,ih)) / gv.vfrac_s;
                                (*gv.c_m).at(il,ih) = ((*gv.c_m).at(il,ih) * gv.vfrac_m 
                                        + (*gv.exchange_si).at(il,ih)) / gv.vfrac_m;
                                (*gv.c_s).at(il,ih) = std::fmax((*gv.c_s).at(il,ih),0.0f);
                                (*gv.c_m).at(il,ih) = std::fmax((*gv.c_m).at(il,ih),0.0f);
                            
                            }else if((*gv.exchange_si).at(il,ih) < 0){
                                //msg = "PROBLEM: negative s->i exchange";
                                //std::cout << msg << std::endl;
                                //print_screen_log(logPULSEfile,&msg);
                            }
                            //if(ih>gv.wetfront_cell){
                            //    (*gv.exchange_si).at(il,ih) = 0.0f;
                            //};
                    }
                };
                //(*gv.c_m) =  ( (*gv.c_m) * gv.vfrac_m_prev + (*gv.exchange_si) * gv.vfrac_s_prev ) / gv.vfrac_m; // / vfrac_m(t);
                //(*gv.c_s) = ( (*gv.c_s) * gv.vfrac_s_prev - (*gv.exchange_si) * gv.vfrac_s_prev ) / gv.vfrac_s;
                
            }
        }
    }

    // Print results                
        if (tcum>=print_next){

            end = std::chrono::system_clock::now();
            elapsed_seconds = end-start;

            outwritestatus = print_results(gv,gp,std::round(print_next),gp.print_step,elapsed_seconds,
                    results_flname);

        if(outwritestatus == true) 
        {
            std::cout << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(gp.print_step/60) << " || Time elapsed (min): " << elapsed_seconds.count()/60 << std::endl;
            (*logPULSEfile) << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(gp.print_step/60) << " || Time elapsed (min): " << std::to_string(elapsed_seconds.count()/60) + "\n";
            print_next += gp.print_step;
            start = std::chrono::system_clock::now();
        } else
        {
            msg = "Problem when saving the results:" + print_next;
            print_screen_log(logPULSEfile,&msg);
            abort();
        }
        }

    }
       
}