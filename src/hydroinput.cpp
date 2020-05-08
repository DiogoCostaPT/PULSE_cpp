
#include "hydroinput.h"

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpPrec(globalvar& gv,double *tcum)
{
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
    
    nqmlt = int((*gv.snowfall_ts).col(0).n_elem) - 1;
    
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