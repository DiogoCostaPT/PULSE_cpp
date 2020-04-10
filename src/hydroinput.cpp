
#include "hydroinput.h"

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpPrec(globalvar& gv,double *tcum)
{
    double prec_t_i = 0.0f,prec_t_i_prev = 0.0f, // time
            prec_i = 0.0f,prec_i_prev=0.0f, // melt rate
            prec_c_i = 0.0f, prec_c_i_prev = 0.0f; // concentration of added snow
    unsigned int nprec = 0;;
    
    nprec = int((*gv.snowfall_ts).col(0).n_elem);
    
    for(int a=0;a<nprec;a++){
        prec_t_i = (*gv.snowfall_ts).at(a,0);
        prec_i = (*gv.snowfall_ts).at(a,1);
        prec_c_i = (*gv.snowfall_ts).at(a,2);
        if(prec_t_i < *tcum){
            prec_t_i_prev = prec_t_i;
            prec_i_prev = prec_i;
            prec_c_i_prev = prec_c_i;
        }else if (prec_t_i == *tcum){
            gv.precip_i =  prec_i;
            break;
        }else if(prec_t_i > *tcum){
            gv.precip_i = prec_i_prev 
                    + (prec_i - prec_i_prev) * ((*tcum - prec_t_i_prev)) 
                    / (prec_t_i - prec_t_i_prev);
            break;
        }
    }
    
    if (gv.precip_i<0){ // acumulation -> add concentration of snow
        if (prec_i<0 && prec_i_prev<0){ // accumulation
            gv.precipc_i = prec_c_i_prev 
                    + (prec_c_i - prec_c_i_prev) * ((*tcum - prec_t_i_prev)) 
                    / (prec_t_i - prec_t_i_prev);
        }else if(prec_i<0 && prec_i_prev>=0){
            gv.precipc_i = prec_c_i;        
        }else if(prec_i>=0 && prec_i_prev<0){
            gv.precipc_i = prec_c_i_prev;        
        }
    }
    
    return;
}

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpQmelt(globalvar& gv,double *tcum)
{
    double qmelt_t_i=0.0f,qmelt_t_i_prev = 0.0f, // time
            qmelt_i=0.0f,qmelt_i_prev=0.0f; // melt rate
    unsigned int nqcmelt = 0;
    
    nqcmelt = int((*gv.qcmel_ts).col(0).n_elem);
    
    for(int a=0;a<nqcmelt;a++){
        qmelt_t_i = (*gv.qcmel_ts).at(a,0);
        qmelt_i = (*gv.qcmel_ts).at(a,1);
        if(qmelt_t_i < *tcum){
            qmelt_t_i_prev = qmelt_t_i;
            qmelt_i_prev = qmelt_i;
        }else if (qmelt_t_i == *tcum){
            gv.qmelt_i =  qmelt_i;
            break;
        }else if(qmelt_t_i > *tcum){
            gv.qmelt_i = qmelt_i_prev 
                    + (qmelt_i - qmelt_i_prev) * ((*tcum - qmelt_t_i_prev)) 
                    / (qmelt_t_i - qmelt_t_i_prev);
            break;
        }
    }
        
    return;
}