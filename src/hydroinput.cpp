
#include "hydroinput.h"

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpMeteo(globalvar& gv,double *tcum)
{
<<<<<<< Updated upstream
    double meteo_t_i = 0.0f,meteo_t_i2 = 0.0f, // time
            temp_calc_i = 0.0f,
            rainfall_calc_i = 0.0f,
            snowfall_calc_i = 0.0f, // melt rate
            precipt_c_calc_i = 0.0f; // concentration of added snow
    unsigned int nprec = 0;;
    
    nprec = int((*gv.meteoall_ts).col(0).n_elem) - 1;
    
    // identification of the time step and linear interpolation for time t
    meteo_t_i2 = 0.0f;
    for(int a=0;a<nprec;a++){
        meteo_t_i = (*gv.meteoall_ts).at(a,0);
        meteo_t_i2 = (*gv.meteoall_ts).at(a+1,0);
        
        if ((*tcum) > meteo_t_i && (*tcum) <= meteo_t_i2) {
            if ((meteo_t_i2 - meteo_t_i) > 0){
                temp_calc_i = (*gv.meteoall_ts).at(a,1);
                rainfall_calc_i = (*gv.meteoall_ts).at(a,2);
                snowfall_calc_i = (*gv.meteoall_ts).at(a,3);
                precipt_c_calc_i = (*gv.meteoall_ts).at(a,4);
                
                gv.tempert_i = temp_calc_i;
                gv.rainfall_i = rainfall_calc_i / (meteo_t_i2 - meteo_t_i); // m/deltatime -> m/sec
                gv.snowfall_i = snowfall_calc_i / (meteo_t_i2 - meteo_t_i); // m/deltatime -> m/sec
                gv.precip_c_i = precipt_c_calc_i;        
             }else
            {
                gv.snowfall_i = 0.0f;
=======
    double prec_t_i = 0.0f,prec_t_i2 = 0.0f, // time
            prec_i = 0.0f, // melt rate
            prec_c_i = 0.0f; // concentration of added snow
    unsigned int nprec = 0;;
    
    nprec = int((*gv.snowfall_ts).col(0).n_elem) - 1;
    
    // identification of the time step and linear interpolation for time t
    prec_t_i2 = 0.0f;
    for(int a=0;a<nprec;a++){
        prec_t_i = (*gv.snowfall_ts).at(a,0);
        prec_t_i2 = (*gv.snowfall_ts).at(a+1,0);
        prec_i = (*gv.snowfall_ts).at(a,1);
        prec_c_i = (*gv.snowfall_ts).at(a,2);

        if ((*tcum) > prec_t_i && (*tcum) <= prec_t_i2) {
            if ((prec_t_i2 - prec_t_i) > 0){
                gv.precip_i = prec_i / (prec_t_i2 - prec_t_i); // m/deltatime -> m/sec
                gv.precipc_i = prec_c_i;        
             }else
            {
                gv.precip_i = 0.0f;
>>>>>>> Stashed changes
            }
        break;
        }
    
    }
}

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpQmelt(globalvar& gv,double *tcum)
{

    double qmlt_t_i = 0.0f,qmlt_t_i2 = 0.0f, // time
            qmlt_i = 0.0f; // melt rate
    unsigned int nqmlt = 0;;
    
<<<<<<< Updated upstream
    nqmlt = int((*gv.meteoall_ts).col(0).n_elem) - 1;
=======
    nqmlt = int((*gv.snowfall_ts).col(0).n_elem) - 1;
>>>>>>> Stashed changes
    
    // identification of the time step and linear interpolation for time t
    qmlt_t_i2 = 0.0f;
    for(int a=0;a<nqmlt;a++){
        qmlt_t_i = (*gv.qcmel_ts).at(a,0);
        qmlt_t_i2 = (*gv.qcmel_ts).at(a+1,0);
        qmlt_i = (*gv.qcmel_ts).at(a,1);

        if ((*tcum) > qmlt_t_i && (*tcum) <= qmlt_t_i2) {
            if ((qmlt_t_i2 - qmlt_t_i) > 0){
                gv.qmelt_i = qmlt_i / (qmlt_t_i2 - qmlt_t_i); // m/deltatime -> m/sec      
             }else
            {
                gv.qmelt_i = 0.0f;
            }
        break;
        }

    }
}