
#include "pulse.h"


/* *****
 * Main PULSE model  
 * **** */
void pulsemodel(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile)
{
    
    // initiation
    double tcum = 0.f,
            //q = 0.f, // melt volume/int
            t = 1.f, 
            deltt = 1.0f, // time step calculated from the CFL condition
            velc_max = 0.0f, // interstitial flow velocity [m s-1]
            D = 0.0f; // dispersion coefficient [m2/s]
    std::string msg;  
    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    bool outwritestatus;
    arma::mat exchange_i;
    double layer_incrmt_l;
    
    unsigned int il,ih,print_next;
      
    tcum = gv.timstart;
    print_next = tcum + gp.print_step;
        
    while (tcum < gp.Tperd)
    {

        t += 1;
                
        findInterpQmelt(gv,&tcum); // if there is increase in SWE, everything will freeze so there will be a stop


        if(gv.q_i==0.0f){ // nothing happens
            tcum++; 
            continue;
        }else if (gv.q_i<0.0f){ // accumulation only 
            upbound_calc(gv,gp,&deltt,logPULSEfile);
            tcum++;
        } else {// melt       
            
            // Estimate interstitial flow velocity 
            (*gv.velc) = gv.q_i / (*gv.vfrac_m); // interstitial flow velocity [m s-1]
            velc_max = arma::max(arma::max(*gv.velc)); 
            (*gv.disp) = gp.aD * (*gv.velc);       // dispersion coefficient [m2/s]
            
            deltt = std::fmin(gp.Courant * gv.snowh / velc_max,gp.print_step);
            
            // limit step so that if there is melt or accumulation it doesn't go more than one cell
            deltt = std::fmin(deltt,gv.snowh*gp.rho_frshsnow_init/(gv.q_i*gp.rho_m));


            tcum = tcum + deltt; 
            
            upbound_calc(gv,gp,&deltt,logPULSEfile);
            
        }

        if (gv.q_i>0.0f){ // if melt
            
            // calculate volume fractions
            vol_fract_calc(gp,gv,&deltt);

            // wetting front
            wetfront_calc(gp,gv,&velc_max,&deltt);

            // check CFC validation
            if((gv.wetfront_cell - gv.wetfront_cell_prev) > 1){
                msg = "CFC condition violation - check code";
                print_screen_log(logPULSEfile,&msg);
                abort();
            }

            // 
            if (arma::min(arma::min(*gv.vfrac_m)) < 1-gp.num_stblty_thrshld_prsity 
                    && arma::min(arma::min(*gv.vfrac_s)) > gp.num_stblty_thrshld_prsity){
              if (gv.wetfront_cell > 5){ // to have sufficient layers for ADE solver

                   crank_nicholson(gv,&deltt); // solve advection and dispersion in the mobile zone

                    // Crank Nicolson to limit the fluxes across boundaries
                   //exchange_i = arma::max(v * deltt * ((*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1) - (*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new))/gv.snowh,(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1));
                   //(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new-1) -= exchange_i; // compute onh advection to the wetting front
                   //(*gv.c_m)(arma::span(0,gv.nl-1),wetfront_cell_new) += exchange_i; // compute onh advection to the wetting front
                }

                (*gv.exchange_si) = (*gv.c_s) * gp.rho_s/gp.rho_m * gv.q_i * deltt / (gv.wetfront_cell+1);

                // limit the flux to the available material
                for(il=0;il<gv.nl ;il++){
                    for(ih=0;ih<gv.nh ;ih++){

                            if ((*gv.exchange_si).at(il,ih) > 0){
                                (*gv.exchange_si).at(il,ih) = std::fmin((*gv.exchange_si).at(il,ih),(*gv.c_s).at(il,ih));
                            }else{
                                msg = "PROBLEM: negative s->i exchange";
                                print_screen_log(logPULSEfile,&msg);
                            }
                            if(ih>gv.wetfront_cell){
                                (*gv.exchange_si).at(il,ih) = 0.0f;
                            };
                    }
                };
                (*gv.c_i) =  ( (*gv.c_i) * gv.vfrac_i_prev + (*gv.exchange_si) % (*gv.vfrac_s_prev) ) / gv.vfrac_i; // / vfrac_m(t);
                (*gv.c_s) = ( (*gv.c_s) % (*gv.vfrac_s_prev) - (*gv.exchange_si) % (*gv.vfrac_s_prev) ) / (*gv.vfrac_s);

                // Exchange with immobile phase (just exchange)
                (*gv.exchange_im)  = deltt * (gp.alphaIE/((*gv.vfrac_m_prev)) % ((*gv.c_i) - (*gv.c_m))); 

                // limit the flux to the available material
                for(il=0;il<gv.nl ;il++){
                    for(ih=0;ih<gv.nh ;ih++){
                         if ((*gv.exchange_im).at(il,ih) > 0){
                            (*gv.exchange_im).at(il,ih) = std::min((*gv.exchange_im).at(il,ih),(*gv.c_i).at(il,ih));
                        }else if((*gv.exchange_im).at(il,ih) < 0){
                             (*gv.exchange_im).at(il,ih) = 0.0f;
                        }
                         //if(ih<gv.upperboundary_cell || ih>gv.wetfront_cell){
                         if(ih>gv.wetfront_cell){
                                (*gv.exchange_im).at(il,ih) = 0.0f;
                            };
                    }
                };
                (*gv.c_m) =  ( (*gv.c_m) % (*gv.vfrac_m_prev) + (*gv.exchange_im) * gv.vfrac_i_prev ) / (*gv.vfrac_m); // / vfrac_m(t);
                (*gv.c_i) = ( (*gv.c_i) * gv.vfrac_i_prev - (*gv.exchange_im) * gv.vfrac_i_prev ) / gv.vfrac_i;

            }
        }

        // Print results                
          if (tcum>=print_next){

             end = std::chrono::system_clock::now();
             elapsed_seconds = end-start;

              outwritestatus = print_results(gv,gp,std::round(print_next),gp.print_step,elapsed_seconds);

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