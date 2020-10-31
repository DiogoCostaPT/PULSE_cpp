
#include "pulse.h"


/* *****
 * Main PULSE model  
 * **** */
bool pulsemodel(globalpar& gp,globalvar& gv,std::ofstream* logPULSEfile,
        std::string* results_flname)
{
    // initiation
    double tcum = 0.f,
            //q = 0.f, // melt volume/int
            t = 0.f,
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

    bool err_flag = false;
        
    tcum = gv.timstart;
    print_next = tcum + gp.print_step;

    while (tcum < gp.Tsim)
    {

        // Snow model: internal or external (SNOWPACK or WARD-snow model)
        if (gp.snowmodel == 0){ // internal    

            // Get precipitation and melt rates            
            findInterpMeteo(gv,&tcum);
            findInterpQmelt(gv,&tcum);

            gv.vfrac_m_prev = gv.vfrac_m;
            gv.vfrac_s_prev = gv.vfrac_s;
            //gv.wetfront_cell_prev = gv.wetfront_cell;

            if (gv.qmelt_t==0.0f && gv.rainfall_t == 0.0f){ // accumulation only 
                deltt = std::fmin(print_next-tcum,
                    gv.v_swe_freshsnow_max/(std::abs(gv.snowfall_t) * gv.snowl));
                velc = 0.0f;
            // watermass_calc_internal(gv,gp,&deltt,&velc,logPULSEfile);
            } else {// melt       
                    // Estimate interstitial flow velocity 
                velc = (gv.qmelt_t + gv.rainfall_t)/ (gv.vfrac_a + gv.vfrac_m); // interstitial flow velocity [m s-1]
                D = gp.aD * velc;       // dispersion coefficient [m2/s]

                (*gv.velc_2d) = (*gv.velc_2d)*0 + velc;
                (*gv.disp_2d) = (*gv.disp_2d)*0 + D;
                
                deltt = std::fmin(print_next-tcum,gp.Courant * gv.snowh / velc);

                // limit step so that if there is melt or accumulation it doesn't go more than one cell
                //deltt_volavail_melt = velc * (*deltt) * max(min((*gv.v_liq)));
                deltt = std::fmin(deltt,gv.v_swe_freshsnow_max/(std::abs(gv.snowfall_t) * gv.snowl));
                deltt = std::fmin(deltt,gv.v_swe_freshsnow_max/(std::abs(gv.qmelt_t) * gv.snowl));
                                
            }

            watermass_calc_internal(gv,gp,&deltt,&velc,logPULSEfile);

            if ((gv.qmelt_t+gv.rainfall_t)>0.0f && gv.nh>0){ // if melt
    
                if (gv.vfrac_m < 1-gp.num_stblty_thrshld_prsity && gv.vfrac_s > gp.num_stblty_thrshld_prsity){  
                    
                    if (gv.wetfront_cell > 3){ // to have sufficient layers for ADE solver

                        if (gp.hydro_solver == 1){
                            crank_nicholson(gv,&deltt,&velc,&D); // solve advection and dispersion in the mobile zone
                        }else if (gp.hydro_solver == 0){
                        //crank_nicholson_hydr2D(gv,&deltt);
                            FtCs_solve_hydr2D(gv,&deltt);

                        }
                        // Limit to available material
                        for(il=0;il<gv.nl ;il++){
                            for(ih=0;ih<gv.nh ;ih++){
                                (*gv.c_m).at(il,ih) = fmax((*gv.c_m).at(il,ih),0.0f);
                            }
                        }
                    }
                }
            }

        } else if (gp.snowmodel == 1){ // external
            
            //  Get time step
            deltt = (*gv.time_ext)(t+1) - (*gv.time_ext)(t); 

            err_flag = watermass_calc_external(gv,gp,&deltt,logPULSEfile, t);
            if (err_flag == true){
                std::string msg = "> ERROR in watermass_calc_external";
                print_screen_log(logPULSEfile,&msg); 
                err_flag = true;
                return err_flag;
            }
            t++;

        } 

        // Ion Exclusion
        if (gv.qmelt_t > 0.0f || gv.rainfall_t > 0.0f){                    
            IonExclusionModel(gp,gv,&deltt);
        }

        // Add one time step
        tcum = tcum + deltt; 

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