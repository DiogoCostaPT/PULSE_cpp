
#include "crank_nicholson_solve_hydr2D.h"

/* *****
* ADE Implicit Solver: Crank-Nicholson Scheme  
* **** */

void crank_nicholson_hydr2D(globalvar& gv,double *deltt)
{
    // calculation - implicit scheme
    int il;    

    // to solveA.x1=B.x0
    int nli = gv.nl;
    int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
    int nt = nli*nhi;
    arma::mat k1 = (*gv.velc_2d)(arma::span(0,nli-1),arma::span(0,nhi-1))*(*deltt)/(4*gv.snowh);       // constants for Crank-Nicholson scheme
    arma::mat k2 = (*gv.disp_2d)(arma::span(0,nli-1),arma::span(0,nhi-1))*(*deltt)/(2*pow(gv.snowh,2));     // constants for Crank-Nicholson scheme
    arma::mat k3 = (*gv.disp_2d)(arma::span(0,nli-1),arma::span(0,nhi-1))*(*deltt)/(2*pow(gv.snowl,2));     // constants for Crank-Nicholson scheme

    k3 *= 0;     // to remove lateral dispersion (onh interested in 1D for now)

    // Creating matrix A to be solved for Crank-Nicholson implicit scheme
    arma::mat a1 = -(k1+k2);
    arma::mat a2 = 1+2*k2+2*k3;
    arma::mat a3 = k1-k2;
    arma::mat a4 =(1-2*k2-2*k3);
    arma::mat a1_r = a1.as_row();
    arma::mat a2_r = a2.as_row();
    arma::mat a3_r = a3.as_row();
    arma::mat a4_r = a4.as_row();
    arma::mat k3_r = k3.as_row();

    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat A(nt,nt,arma::fill::zeros);
    A +=  arma::diagmat(a3_r(0,arma::span(nli,nt-1)),nli)
          +arma::diagmat(a1_r(0,arma::span(0,nt-nli-1)),-nli)
          +arma::diagmat(a2_r(0,arma::span(0,nt-1)))
          -arma::diagmat(k3_r(0,arma::span(1,nt-2)),1)
          //-arma::diagmat(k3_r(0,arma::span(2,nt-1)),-1)
        ;

    //  I think it's creting the higher values in the margins
    //for(il=0;il<(nt/nli)-1 ;il++){ // inner cells - left and right marigns
    //    A(il*nli+1,il*nli)=0;
    //    A(il*nli+1,il*nli+1)=A(il*nli+1,il*nli+1)-k3(il*nli+1,il*nli+1);
    //    A(il*nli+nli,il*nli+nli+1)=0;
    //    A(il*nli+nli,il*nli+nli)=A(il*nli+nli,il*nli+nli)-k3(il*nli+1,il*nli+1);
    //}

    // 2) first row (y=1)
    A(arma::span(0,nli-1),arma::span(0,nli-1)) = // check if shoudn't be +=
        arma::diagmat(a1_r); //diagonal
    A(0,0)=a1(0,0)+a2(0,0)-k3(0,0); // left corner
    A(nli-1,nli-1)=a2(nli-1,0)+a1(nli-1,0)-k3(nli-1,0); // left corner

    
    // 3) last row (y=ny) 
    A(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1)) = 
        arma::diagmat(a2_r + a3_r); //diagonal 
    A(nt-1,nt-1)=a2(0,nhi-1)+a3(0,nhi-1)-k3(0,nhi-1); // right corner
    A(nt-nli,nt-nli)=a3(nli-1,nhi-1)+a2(nli-1,nhi-1)-k3(nli-1,nhi-1); // left corner

    // Creating B


    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat B(nt,nt,arma::fill::zeros);
    B = -arma::diagmat(a3_r(0,arma::span(0,nt-nli-1)),nli)
            -arma::diagmat(a1_r(0,arma::span(0,nt-nli-1)),-nli)
            +arma::diagmat(a4_r(0,arma::span(0,nt-1)))
            +arma::diagmat(k3_r(0,arma::span(0,nt-2)),1);
            +arma::diagmat(k3_r(0,arma::span(0,nt-2)),-1);


    //for(il=0;il<(nt/nli)-1 ;il++){         // inner cells - left and right marigns
    //   B(il*nli+1,il*nli)=0;
    //    B(il*nli+1,il*nli+1)=A(il*nli+1,il*nli+1)-k3(il*nli+1,il*nli+1);
    //    B(il*nli+nli,il*nli+nli+1)=0;
    //    B(il*nli+nli,il*nli+nli)=A(il*nli+nli,il*nli+nli)-k3(il*nli+1,il*nli+1);
    //}

    // 2) first row (y=1)
    B(arma::span(0,nli-1),arma::span(0,nli-1)) = 
        arma::diagmat(-a1_r); //diagonal
    B(0,0)=-a1(0,0)+a4(0,0)+k3(0,0); // left corner
    B(nli-1,nli-1)=a4(nli-1,0)-a1(nli-1,0)+k3(nli-1,0); // left corner


    // 3) last row (y=ny)
    B(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1)) = 
        arma::diagmat(a4_r - a3_r); //diagonal   
    B(nt-1,nt-1)=a4(0,nhi-1)-a3(0,nhi-1)+k3(0,nhi-1); // right corner
    B(nt-nli,nt-nli)=-a3(nli-1,nhi-1)+a4(nli-1,nhi-1)+k3(nli-1,nhi-1); // left corner

    arma::mat c0 = (*gv.c_m)(arma::span(0,nli-1),arma::span(0,(gv.wetfront_cell)-1));
    arma::mat c1 = arma::reshape(c0,1,nt);  //c1=reshape(c_m_prev(2:end-1,2:end-1),1,[]);
    arma::vec b=B*trans(c1);    // calculation of [b]
    arma::mat c2=arma::solve(A,b);     // calculation of c
    (*gv.c_m)(arma::span(0,nli-1),arma::span(0,gv.wetfront_cell-1)) = arma::reshape(c2,nli,nhi);

}