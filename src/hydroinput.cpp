
#include "hydroinput.h"

/* *****
 * Determine snowmelt rate and, in the case of accumulation, determine the concentration of the added snow layer 
 * **** */
void findInterpQmelt(globalvar& gv,double *tcum)
{
    double qcmelt_t_i,qcmelt_t_i_prev = 0.0f, // time
            qcmelt_i,qcmelt_i_prev=0.0f, // melt rate
            qcmelt_c_i = 0.0f, qcmelt_c_i_prev = 0.0f; // concentration of added snow
    unsigned a,nqcmelt;
    
    nqcmelt = int((*gv.qcmel_ts).col(0).n_elem);
    
    for(a=0;a<nqcmelt;a++){
        qcmelt_t_i = (*gv.qcmel_ts).at(a,0);
        qcmelt_i = (*gv.qcmel_ts).at(a,1);
        qcmelt_c_i = (*gv.qcmel_ts).at(a,2);
        if(qcmelt_t_i < *tcum){
            qcmelt_t_i_prev = qcmelt_t_i;
            qcmelt_i_prev = qcmelt_i;
            qcmelt_c_i_prev = qcmelt_c_i;
        }else if (qcmelt_t_i == *tcum){
            gv.q_i =  qcmelt_i;
            break;
        }else if(qcmelt_t_i > *tcum){
            gv.q_i = qcmelt_i_prev 
                    + (qcmelt_i - qcmelt_i_prev) * ((*tcum - qcmelt_t_i_prev)) 
                    / (qcmelt_t_i - qcmelt_t_i_prev);
            break;
        }
    }
    
    if (gv.q_i<0){ // acumulation -> add concentration of snow
        if (qcmelt_i<0 && qcmelt_i_prev<0){ // accumulation
            gv.qc_i = qcmelt_c_i_prev 
                    + (qcmelt_c_i - qcmelt_c_i_prev) * ((*tcum - qcmelt_t_i_prev)) 
                    / (qcmelt_t_i - qcmelt_t_i_prev);
        }else if(qcmelt_i<0 && qcmelt_i_prev>=0){
            gv.qc_i = qcmelt_c_i;        
        }else if(qcmelt_i>=0 && qcmelt_i_prev<0){
            gv.qc_i = qcmelt_c_i_prev;        
        }
    }
    
    return;
}