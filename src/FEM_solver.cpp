
#include "FEM_solver.h"

/* *****
* ADE Implicit Solver: Crank-Nicholson Scheme  
* **** */

void FtCs_solve_hydr2D(globalvar& gv,double *deltt)
{

int nli = gv.nl;
int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
int nt = nli*nhi;
int ih,il;

arma::mat c_t1 = (*gv.c_m)(arma::span(0,nli-1),arma::span(0,(gv.wetfront_cell)-1));
arma::mat c_t2 = c_t1;

arma::mat dcdy(nli,gv.wetfront_cell,arma::fill::zeros); 
arma::mat dc2dy2(nli,gv.wetfront_cell,arma::fill::zeros); 
arma::mat dc2dx2(nli,gv.wetfront_cell,arma::fill::zeros); 
arma::mat dcdt(nli,gv.wetfront_cell,arma::fill::zeros); 

for (ih=1;ih<gv.wetfront_cell-1;ih++){
    for (il=1;il<nli-1;il++){
        
        if ((*gv.velc_2d)(il,ih-1) > 0){
            //dtcm = ((c1(il,ih-1) * gv.vfrac_m) * (*gv.velc_2d)(il,ih-1) * (*deltt))/gv.vfrac_m;
            //dtcm = std::max(dtcm,c1(il,ih-1));
            //c1(il,ih-1) -= dtcm;
            //c1(il,ih) += dtcm;

            dcdy(il,ih)=(c_t1(il,ih)-c_t1(il,ih-1))/gv.snowh;
            dc2dy2(il,ih)=((c_t1(il,ih+1)-c_t1(il,ih))/gv.snowh-(c_t1(il,ih)-c_t1(il,ih-1))/gv.snowh)/(gv.snowh);         
            dc2dx2(il,ih)=((c_t1(il+1,ih)-c_t1(il,ih))/gv.snowl-(c_t1(il,ih)-c_t1(il-1,ih))/gv.snowl)/(gv.snowl);
            dcdt(il,ih)=-(*gv.velc_2d)(il,ih)*dcdy(il,ih)+(*gv.disp_2d)(il,ih)*(dc2dx2(il,ih)+dc2dy2(il,ih));

            c_t2(il,ih)=c_t2(il,ih)+dcdt(il,ih)*(*deltt);

        }
    }
}

// Von Neumann condition: dc/dx=0 in the river margins
c_t2(0,arma::span(0,nhi-1))=c_t2(1,arma::span(0,nhi-1));  
c_t2(nli-1,arma::span(0,nhi-1))=c_t2(nli-2,arma::span(0,nhi-1));   
    
(*gv.c_m)(arma::span(0,nli-1),arma::span(0,gv.wetfront_cell-1)) =  c_t2;

};