// Copyright 2021: Diogo Costa

// This program, PULSE_cpp, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "crank_nicholson_solve.h"

/* *****
 * ADE Implicit Solver: Crank-Nicholson Scheme  
 * **** */
void crank_nicholson(globalvar& gv,double *deltt,double *v,double *D)
{
    // calculation - implicit scheme
    int il;    

    // to solveA.x1=B.x0
    int nli = gv.nl;
    int nhi = gv.wetfront_cell;//-(gv.upperboundary_cell);                   // the boundaries are knowns, so don't need to be included in matrix A
    int nt = nli*nhi;
    double k1 = (*v)*(*deltt)/(4*gv.snowh);       // constants for Crank-Nicholson scheme
    double k2 = (*D)*(*deltt)/(2*pow(gv.snowh,2));     // constants for Crank-Nicholson scheme
    double k3 = (*D)*(*deltt)/(2*pow(gv.snowl,2));     // constants for Crank-Nicholson scheme

    k3=0;     // to remove lateral dispersion (onh interested in 1D for now)

    // Creating matril A to be solved for Crank-Nicholson implicit scheme
    double a1=-(k1+k2);
    double a2=(1+2*k2+2*k3);
    double a3=k1-k2;

    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat A(nt,nt,arma::fill::zeros);
    A = a3*(arma::diagmat(arma::ones(1,nt-nli),nli))
                    +a1*arma::diagmat(arma::ones(1,nt-nli),-(nli))
                    +a2*arma::diagmat(arma::ones(1,nt))
                    -k3*arma::diagmat(arma::ones(1,nt-1),1)
                    -k3*arma::diagmat(arma::ones(1,nt-1),-1);


    for(il=0;il<(nt/nli)-1 ;il++){ // inner cells - left and right marigns
        A(il*nli+1,il*nli)=0;
        A(il*nli+1,il*nli+1)=A(il*nli+1,il*nli+1)-k3;
        A(il*nli+nli,il*nli+nli+1)=0;
        A(il*nli+nli,il*nli+nli)=A(il*nli+nli,il*nli+nli)-k3;
    }

    // 2) first row (y=1)
    A(arma::span(0,nli-1),arma::span(0,nli-1)) = A(arma::span(0,nli-1),arma::span(0,nli-1)) 
                                                + a1*arma::diagmat(arma::ones(1,nli)); //diagonal
    A(0,0)=a1+a2-k3; // left corner
    A(nli-1,nli-1)=a2+a1-k3; // left corner

    
    // 3) last row (y=ny)
    A(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1)) = A(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1))
                            + (a2+a3)*arma::diagmat(arma::ones(1,nli)); // diagonal !!!! CHECK IF IT SHOULDN'T BE "+=" (LINE ABOVE) INSTEAD OF "="
    A(nt-1,nt-1)=a2+a3-k3; // right corner
    A(nt-nli,nt-nli)=a3+a2-k3; // left corner

    // Creating B
    double a4=(1-2*k2-2*k3);

    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat B(nt,nt,arma::fill::zeros);
    B=-a3*arma::diagmat(arma::ones(1,nt-nli),nli)
            -a1*arma::diagmat(arma::ones(1,nt-nli),-(nli))
            +a4*arma::diagmat(arma::ones(1,nt))
            +k3*arma::diagmat(arma::ones(1,nt-1),1)
            +k3*arma::diagmat(arma::ones(1,nt-1),-1);


    for(il=0;il<(nt/nli)-1 ;il++){         // inner cells - left and right marigns
        B(il*nli+1,il*nli)=0;
        B(il*nli+1,il*nli+1)=B(il*nli+1,il*nli+1)-k3;
        B(il*nli+nli,il*nli+nli+1)=0;
        B(il*nli+nli,il*nli+nli)=B(il*nli+nli,il*nli+nli)-k3;
    }

    // 2) first row (y=1)
    B(arma::span(0,nli-1),arma::span(0,nli-1)) = B(arma::span(0,nli-1),arma::span(0,nli-1)) 
                                                + (-a1)*arma::diagmat(arma::ones(1,nli)); //diagonal
    B(0,0)=-a1+a4+k3; // left corner
    B(nli-1,nli-1)=a4-a1+k3; // left corner

    // 3) last row (y=ny)
    B(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1)) = B(arma::span(nt-nli,nt-1),arma::span(nt-nli,nt-1))
                        +(a4-a3)*arma::diagmat(arma::ones(1,nli)); // diagonal !!!! CHECK IF IT SHOULDN'T BE "+=" (LINE ABOVE) INSTEAD OF "="
    B(nt-1,nt-1)=a4-a3+k3; // right corner
    B(nt-nli,nt-nli)=-a3+a4+k3; // left corner

    arma::mat c0 = (*gv.c_m)(arma::span(0,nli-1),arma::span(0,(gv.wetfront_cell)-1));
    arma::mat c1 = arma::reshape(c0,1,nt);  //c1=reshape(c_m_prev(2:end-1,2:end-1),1,[]);
    arma::vec b=B*trans(c1);    // calculation of [b]
    arma::mat c2=arma::solve(A,b);     // calculation of c
    (*gv.c_m)(arma::span(0,nli-1),arma::span(0,gv.wetfront_cell-1)) = arma::reshape(c2,nli,nhi);


}
