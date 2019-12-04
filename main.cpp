
// Copyright 2017-2019: Diogo Costa

// This program, PULSE, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) an_col later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include<iostream>
#include<fstream>
#include<math.h>
#include<armadillo>
#include<string>
#include<memory> 
#include <chrono>
#include <ctime>  

#include <vector> 
#include <dirent.h>
#include <sys/types.h>


class globalpar
{
public:
    
    double Courant=0.8,aD,alpha,rho_s,rho_m,
           wetfront_z,num_stblty_thrshld_prsity = 1E-6,ionexcl_fact,Tperd;
    
    int flag_sens,run_id,s,print_step;
    //std::ofstream logPULSEfile;
    
};


class globalvar
{
public:
  globalvar() {

  }
  globalvar(size_t nh, size_t nl) {
    this->nh = nh;
    this->nl = nl;
   
    c_m = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nh,nl));
    c_i = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nh,nl));
    c_s = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nh,nl));
    exchange_si = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nh,nl));
    exchange_im = std::unique_ptr<arma::Mat<double>>( new  arma::mat(nh,nl));
    
    qmelt = std::unique_ptr<arma::Mat<double>>( new  arma::mat(2000,2));

  }
    
    size_t nh,nl,
        snowH, // snowpack depth
        snowL, // snowpack horizontal lenght
        snowl, // grid h lenght
        snowh; // grtid l lenght
  
    std::unique_ptr<arma::Mat<double>> c_m,c_i,c_s,qmelt,exchange_si,exchange_im;
    double porosity_m=0.008,
            porosity_i=0.001,
            porosity_s=1 - porosity_m - porosity_i,
            porosity_m_prev=porosity_m,
            porosity_i_prev=porosity_i,
            porosity_s_prev=porosity_s;
                 
};



void print_screen_log(std::ofstream& logPULSEfile,std::string& msg)
{
    std::cout << msg << std::endl;
    logPULSEfile << msg + "\n";  
    
}

void checkmesh(int* H_local,int* L_local,int* h_layer,int* l_layer,int* nl,int* nh,std::ofstream&logPULSEfile)
{
    float a,b;
    std::string msg;
    
    a = (*H_local)/(*h_layer);
    b = (*L_local)/(*l_layer);
    
    if(floor(a)==ceil(a) && floor(b)==ceil(b)){      
        (*nh) = a;
        (*nl) = b;
        msg = "Snowpack mesh: created";
        print_screen_log(logPULSEfile,msg);
    }else{
        if (floor(a)==ceil(a)){
        msg = "Snowpack mesh: H" + std::to_string ((*H_local)) + 
                "(snowpack depth) is not divisible by h_layer" + std::to_string((*h_layer)) +
                "(grid thickness) -> change the simulation setting in file 'simset.pulse'";
        print_screen_log(logPULSEfile,msg);
        }
        if (floor(b)==ceil(b)){
        msg = "Snowpack mesh: L" + std::to_string ((*L_local)) + 
                "(snowpack horizontal length) is not divisible by l_layer" + std::to_string((*l_layer)) +
                "(grid horizontal length) -> change the simulation setting in file 'simset.pulse'";
        print_screen_log(logPULSEfile,msg);
        }
        abort();
    }
    
}
    

void read_simset(globalpar& gp,std::string& sim_purp, int *H_local,int *L_local, int *h_layer,int *l_layer, std::string& qmelt_file,std::ofstream& logPULSEfile)
{
    // read_modset(ds,print_step,ks_input,zbinc,ntim_days)
    
    std::string str, modset_flname, msg;
    modset_flname = "simset.pulse";
    
    std::ifstream file(modset_flname);
    
    int i = 0;
    while (std::getline(file, str)) 
    {
        i += 1;
        if(i==1){sim_purp = str;};
        if(i==2){(*H_local) = std::stoi(str);};
        if(i==3){(*L_local) = std::stoi(str);};
        if(i==4){(*h_layer) = std::stoi(str);};
        if(i==5){(*l_layer) = std::stoi(str);};
        if(i==6){qmelt_file = str;};  
        if(i==7){gp.print_step = std::stoi(str);};  
    }
    file.close();
    
    if(i==7){
        msg = "Successful loading the file: " + modset_flname;
    } else{
        msg = "PROBLEM loading the file: " + modset_flname;
    } 
    print_screen_log(logPULSEfile,msg);   
    
}    


void read_qmelt(globalvar& gv,std::string& qmelt_file,std::ofstream& logPULSEfile)
{
    unsigned int a; 
    double tmelts,vmelt, tmelts_bef = 0.0f;
    gv.qmelt  = 0;
    std::string msg;
    
    arma::mat filedataQ; 
    bool flstatusQ =  filedataQ.load(qmelt_file,arma::csv_ascii);
    if(flstatusQ == true) {
        for(a=0;a<filedataQ.col(1).n_elem;a++){
            tmelts = filedataQ(a,0);  // t melt seconds
            vmelt = filedataQ(a,1);  // value of melt
            (*gv.qmelt).at(a,0) = tmelts;  
            (*gv.qmelt).at(a,1) = vmelt;
            //gv.qmelt += vmelt /(1000.*3600.*24.) * (tmelts - tmelts_bef); 
            tmelts_bef = tmelts;
        }
       msg = "Successful loading the file: " + qmelt_file;
    } else{
        msg = "PROBLEM loading the file: " + qmelt_file;       
    } 
    print_screen_log(logPULSEfile,msg);   
    
    return;
  
}

void calc_porosity(globalpar& gp,globalvar& gv,double *q, double *deltt)
{
    
    double dporosity_s_dt = (*q) / arma::accu((*gv.qmelt));
    double dporosity_i_dt = dporosity_s_dt;

    gv.porosity_m_prev = gv.porosity_m;
    gv.porosity_i_prev = gv.porosity_i;
    gv.porosity_s_prev = gv.porosity_s;
    
    gv.porosity_s = std::fmax(gv.porosity_s - dporosity_s_dt * (*deltt), 0.f);
    if (gv.porosity_s != 0){
        gv.porosity_i = std::fmax(gv.porosity_i + dporosity_s_dt * (*deltt) - dporosity_i_dt * (*deltt), 0.f); 
    }else{
        gv.porosity_i = std::fmax(gv.porosity_i - dporosity_i_dt * (*deltt) , 0.f);     
    };
    gv.porosity_m = std::fmin(gv.porosity_m + dporosity_i_dt * (*deltt) , 1.f);  
       
}

void wettingfront_cell_location(globalpar& gp,globalvar& gv,double *v, double *deltt, int *wetfront_cell_new,int *wetfront_cell_prev)
{
    int nh_l = gv.nh;
    *wetfront_cell_prev = *wetfront_cell_new; // wetting fron cell in the previous time step
    gp.wetfront_z = std::fmax(gp.wetfront_z - (*v) * (*deltt),0.f);
    int tmp = std::round(nh_l-gp.wetfront_z/gv.snowh+1);
    *wetfront_cell_new = std::min(tmp,nh_l); // finding the cell when the wetting front is located
    //wetfront_cell = (*wetfront_cell_new);
    
}


// Crank-Nicholson Scheme (implicit)
void Crank_Nicholson(globalvar& gv,int *wetfront_cell_prev, int *wetfront_cell_new, double *deltt,double *v,double *D)
{
    // calculation - implicit scheme
    unsigned int i,j;    

    // to solve A.x1=B.x0
    int nxi = *wetfront_cell_new-2;
    int nyi = gv.nh-2;                  // the boundaries are knowns, so don't need to be included in matrix A
    int nt = nxi*nyi;
    double k1 = (*v)*(*deltt)/(4*gv.snowh);       // constants for Crank-Nicholson scheme
    double k2 = (*D)*(*deltt)/(2*pow(gv.snowh,2));     // constants for Crank-Nicholson scheme
    double k3 = (*D)*(*deltt)/(2*pow(gv.snowl,2));     // constants for Crank-Nicholson scheme

    k3=0;     // to remove lateral dispersion (onh interested in 1D for now)

    // Creating matrix A to be solved for Crank-Nicholson implicit scheme
    double a1=-(k1+k2);
    double a2=(1+2*k2+2*k3);
    double a3=k1-k2;
    //A=zeros(nt,nt);

    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat A = a3*(arma::diagmat(arma::ones(1,nt-nxi),nxi))+
                    a1*arma::diagmat(arma::ones(1,nt-nxi),-(nxi))+
                    a2*arma::diagmat(arma::ones(1,nt))
                    -k3*arma::diagmat(arma::ones(1,nt-1),1)-
                    k3*arma::diagmat(arma::ones(1,nt-1),-1);


    for(i=1;i<(nt/nxi)-2 ;i++){ // inner cells - left and right marigns
        A(i*nxi+1,i*nxi)=0;
        A(i*nxi+1,i*nxi+1)=A(i*nxi+1,i*nxi+1)-k3;
        A(i*nxi+nxi,i*nxi+nxi+1)=0;
        A(i*nxi+nxi,i*nxi+nxi)=A(i*nxi+nxi,i*nxi+nxi)-k3;
    }

    // 2) first row (y=1)
    A(arma::span(1,nxi),arma::span(1,nxi)) += a1*arma::diagmat(arma::ones(1,nxi)); //diagonal
    A(1,1)=a1+a2-k3; // left corner
    A(nxi,nxi)=a2+a1-k3; // left corner
    A(nxi,nxi+1)=0; // left corner

    
    // 3) last row (y=ny)
    A(arma::span(nt-nxi+1,nt),arma::span(nt-nxi+1,nt))=(a2+a3)*arma::diagmat(arma::ones(1,nxi)); // diagonal
    A(nt,nt)=a2+a3-k3; // right corner
    A(nt-nxi+1,nt-nxi+1)=a3+a2-k3; // left corner
    A(nt-nxi+1,nt-nxi)=0; // left corner


    // Creating B
    double a4=(1-2*k2-2*k3);
    //B=zeros(nt,nt);

    // 1) all domain (diagonals: -9,-1,0,1,9)
    arma::mat B=-a3*arma::diagmat(arma::ones(1,nt-nxi),nxi)-
            a1*arma::diagmat(arma::ones(1,nt-nxi),-(nxi))+
            a4*arma::diagmat(arma::ones(1,nt))+
            k3*arma::diagmat(arma::ones(1,nt-1),1)+
            k3*arma::diagmat(arma::ones(1,nt-1),-1);


    for(i=1;i<(nt/nxi)-2 ;i++){         // inner cells - left and right marigns
        B(i*nxi+1,i*nxi)=0;
        B(i*nxi+1,i*nxi+1)=B(i*nxi+1,i*nxi+1)+k3;
        B(i*nxi+nxi,i*nxi+nxi+1)=0;
        B(i*nxi+nxi,i*nxi+nxi)=B(i*nxi+nxi,i*nxi+nxi)+k3;
    }

    // 2) first row (y=1)
      
    B(arma::span(1,nxi),arma::span(1,nxi)) += (-a1)*arma::diagmat(arma::ones(1,nxi)); //diagonal
    B(1,1)=-a1+a4+k3; // left corner
    B(nxi,nxi)=a4-a1+k3; // left corner
    B(nxi,nxi+1)=0; // left corner

    // 3) last row (y=ny)
    B(arma::span(nt-nxi+1,nt),arma::span(nt-nxi+1,nt))=(a4-a3)*arma::diagmat(arma::ones(1,nxi)); // diagonal
    B(nt,nt)=a4-a3+k3; // right corner
    B(nt-nxi+1,nt-nxi+1)=-a3+a4+k3; // left corner
    B(nt-nxi+1,nt-nxi)=0; // left corner

    //c_m_new = c_m_prev;
    
    int nxyi_act = (nxi-2) * (nyi-2);
    arma::mat c1 = arma::reshape((*gv.c_m)(arma::span(2,nxi-1),arma::span(2,nyi-1)),1,nxyi_act);  //c1=reshape(c_m_prev(2:end-1,2:end-1),1,[]);
    arma::mat b=B*trans(c1);    // calculation of [b]
    arma::mat c2=arma::solve(A,b);     // calculation of c
    (*gv.c_m)(arma::span(2,nxi-1),arma::span(2,nyi-1)) = arma::reshape(c2,nxi,nyi);


    for(i=1;i<nxi ;i++){ // to avoid back flow due to roundoff errors
        for(j=1;j<nyi ;j++){ 
             if ((*gv.c_m)(i,j)<0){
                 (*gv.c_m)(i,j)=0;
            }
        }
    }
}


bool print_results(globalvar& gv,globalpar& gp, int print_tag, unsigned int print_step, std::chrono::duration<double> elapsed_seconds)
{

    unsigned int il,ih;
    int a = 0;
    double ux;
    
    std::string tprint = "Results/" + std::to_string(print_tag); 
    std::string filext(".txt");
    tprint += filext;

    arma::mat filedataR(gv.nl*gv.nh,9); 
    
    for(ih=1;ih<=gv.nh;ih++)
    {
        for(il=1;il<=gv.nl;il++)
        {
            filedataR(a,0) = ih * gv.snowh;  
            filedataR(a,1) = il * gv.snowl;  
            filedataR(a,2) = (*gv.c_m).at(ih,il); 
            filedataR(a,3) = (*gv.c_i).at(ih,il); 
            filedataR(a,4) = (*gv.c_s).at(ih,il); 
            filedataR(a,5) = (gv.porosity_m); 
            filedataR(a,6) = (gv.porosity_s); 
            filedataR(a,7) = (*gv.exchange_si).at(ih,il); 
            filedataR(a,8) = (*gv.exchange_im).at(ih,il); 
            a = a + 1;
        }
    }
   
    arma::mat filedata(std::max(0,a-1),8); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,8));
    
    bool outwritestatus =  filedata.save(tprint,arma::csv_ascii);
    return outwritestatus;
}

void PULSEmodel(globalpar& gp,globalvar& gv,std::ofstream& logPULSEfile)
{
    
    // initiation
    double tcum = 0.f,
            q = 0.f, // melt volume/int
            t = 1.f, 
            deltt = 1.0f, // time step calculated from the CFL condition
            v = 0.f, // interstitial flow velocity [m s-1]
            D = 0.f; // dispersion coefficient [m2/s]
    int upperboundary_cell = 1,
        wetfront_cell_new,wetfront_cell_prev,upperboundary_cell_new,upperboundary_cell_prev,
        flagt = 1, // for saving results
        snowh_min,tmp_int; // min vertical grid size to comply with the Peclet condition
    std::string msg;  
    double upperboundary_z = gp.wetfront_z,exchange_i;
    double Peclet,Peclet_max = 0;
    std::chrono::duration<double> elapsed_seconds;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    bool outwritestatus;
    
    gp.wetfront_z = gv.nh*gv.snowh; // starts at the top
    unsigned int ix,iy, nh_l = gv.nh, nl_l = gv.nl,print_next;
    
    print_next = gp.print_step;
    while (tcum < gp.Tperd)
    {

        t += 1;

        q = std::fmax((*gv.qmelt).at(floor(tcum+1)),0); // if there is increse in SWE, everything will freeze so there will be a stop

        // Estimate interstitial flow velocity 
        v = q / (gv.porosity_m); // interstitial flow velocity [m s-1]
        D = gp.aD * v;       // dispersion coefficient [m2/s]

        // Dynamic time interval to comply with Courant Condition
        if (q != 0){
            deltt = std::fmin(gp.Courant * gv.snowh / v,1);
        }else{
            deltt = 1;
        }

        // Calculate the Peclet number
        if (v > 0){
            Peclet = (v * gv.snowh)/D;
            if (Peclet > 2 && Peclet > Peclet_max){
                snowh_min = 2 * D / v;
                msg = "Peclet number > 2. Delta y needs to be equal or smaller than " + std::to_string(ceil(snowh_min));
                print_screen_log(logPULSEfile,msg); 
            }
            Peclet_max = Peclet;
        }
        tcum = tcum + deltt; 

        // Calculate porosity for next time step
        calc_porosity(gp,gv,&q,&deltt);
        
        // limiting the flux to the wetting front (uses intersticial velocity to determine the wetting front)
        wettingfront_cell_location(gp,gv,&v,&deltt,&wetfront_cell_new,&wetfront_cell_prev);
        
        // Getting the active layers
        //-// arma::mat c_m_prev_i = (*gv.c_m).rows(1,wetfront_cell_prev); // moved inside the solver now in cpp
        //-// arma::mat c_i_prev_i = (*gv.c_i).rows(1,wetfront_cell_prev); // moved inside the solver now in cpp
        //-// arma::mat c_s_prev_i = (*gv.c_s).rows(1,wetfront_cell_prev); // moved inside the solver now in cpp
        //c_m_new_i=c_m(:,1:wetfront_cell,t);
        //c_i_new_i=c_i(:,1:wetfront_cell,t);
        //c_s_new_i=c_s(:,1:wetfront_cell,t);

        // Melting velocity: last cell
        //upperboundary_cell_prev = upperboundary_cell; // wetting fron cell in the previous time step
        //upperboundary_z =  std::fmax(upperboundary_z - q * deltt,0);
        //tmp_int = std::round(nh_l-upperboundary_z/gv.snowh+1);
        //upperboundary_cell_new = std::min(tmp_int,nh_l); // in what cell is the wetting front
        //upperboundary_cell = [upperboundary_cell, upperboundary_cell_new];
        
        
        // add a new cell if the wetting front moves to the next cell
        if(wetfront_cell_new - wetfront_cell_prev > 1){
            msg = "Courant condition violation - check code";
            print_screen_log(logPULSEfile,msg); 

            break;
        //-// }else{ // moved inside the solver now in cpp
        //-//     if (wetfront_cell_new == wetfront_cell_prev + 1){ // moved inside the solver now in cpp
             //-// c_m_prev_i =  arma::join_vert(c_m_prev_i,arma::zeros(nl_l,1));  // moved inside the solver now in cpp
             //c_i_prev_i = [c_i_prev_i(:,wetfront_cell_prev),zeros(nl,1)]; 
             //c_s_prev_i = [c_s_prev_i(:,wetfront_cell_prev),ones(nl,1)*c0_s];
             //-// arma::mat c_m_new_i = c_m_prev_i; // moved inside the solver now in cpp

             //arma::mat c_i_new_i = [c_i_prev_i(:,1:wetfront_cell_prev),zeros(nl,1)]; 
             //arma::mat c_s_new_i = [c_s_prev_i(:,1:wetfront_cell_prev),c_s(:,wetfront_cell_prev+1)];

             //-// arma::mat c_i_new_i = arma::join_vert(c_i_prev_i,arma::zeros(nl_l,1));  // moved inside the solver now in cpp
             //-// arma::mat c_s_new_i = arma::join_vert(c_s_prev_i,(*gv.c_s).row(wetfront_cell_prev+1));  // moved inside the solver now in cpp
           
        //-//     }else{ // moved inside the solver now in cpp
        //-//      arma::mat c_s_new_i = c_s_prev_i;
        //-//      arma::mat c_i_new_i = c_i_prev_i;
        //-//      arma::mat c_m_new_i = c_m_prev_i;
        //-//  } 
        }

        if (gv.porosity_m < 1-gp.num_stblty_thrshld_prsity && gv.porosity_i > gp.num_stblty_thrshld_prsity && gv.porosity_s > gp.num_stblty_thrshld_prsity){
            if (wetfront_cell_new > 5){
              //int nh_Crank = wetfront_cell_new - 1;        
              //_m_prev_ii = c_m_prev_i(:,1:nh_Crank);
              //c_m_prev_ii = c_m_prev_i(:,1:wetfront_cell_new);

              // c_m_comp = Crank_Nicholson(t,nl,wetfront_cell_new,c_m_prev_ii,snowl,snowh,deltt,v(t),D(t)); % solve advection and dispersion in the mobile zone
              Crank_Nicholson(gv,&wetfront_cell_prev, &wetfront_cell_new,&deltt,&v,&D); // solve advection and dispersion in the mobile zone

              // Computing the wet front - it needs to be calculated seperatly from
              // Crank Nicolson to limit the fluxes across boundaries

             (*gv.c_m)(arma::span(0,gv.nh),wetfront_cell_new) = (*gv.c_m)(arma::span(0,gv.nh),wetfront_cell_new) - v * deltt * ((*gv.c_m)(arma::span(0,gv.nh),wetfront_cell_new) - (*gv.c_m)(arma::span(0,gv.nh),wetfront_cell_new-1))/gv.snowh ; // compute onh advection to the wetting front

             //c_m_new_i(:,wetfront_cell_new) = c_m_new_i_wetfront;

             // fluxback = sum(c_m(:,wetfront_cell+1:end,t)');
             //c_m(:,wetfront_cell-1,t) = c_m(:,wetfront_cell-1,t) + fluxback';

            //-//}else{
                // do nothing, onh add the sources
            //-//    c_m_new_i =  c_m_prev_i;
            }

            // c_i_new_i = c_i_prev_i;
            // c_s_new_i = c_s_prev_i;

            // Effect of freeze-thaw
            // exchange = (c_s_new_i -  c_m_new_i) * rho_s/rho_m * q; %dporosity_s_dt;
            // exchange = c_s_new_i * rho_s/rho_m ; %dporosity_s_dt;

            (*gv.exchange_si) = (*gv.c_s)* gp.rho_s/gp.rho_m * q * deltt / wetfront_cell_new * gp.ionexcl_fact;

            // limit the flux to the available material
            for(ix=1;ix<gv.nl ;ix++){
                for(iy=1;iy<gv.nh ;iy++){
                    if ((*gv.exchange_si)(ix,iy) > 0){
                        (*gv.exchange_si)(ix,iy) = std::min((*gv.exchange_si)(ix,iy),(*gv.c_s)(ix,iy));
                    }else if((*gv.exchange_si)(ix,iy) < 0){
                        msg = "PROBLEM: negative exchange";
                        print_screen_log(logPULSEfile,msg);
                        //-// (*gv.exchange)(ix,iy)  = - std::min(std::abs(exchange),abs(c_m_new_i)); 
                    }
                }
            };

            // c_m(:,wetfront_cell,t) = c_m(:,wetfront_cell,t) + u * c_s(:,wetfront_cell,t) / rho_m;
            // c_i_new_i =  c_i_new_i + exchange * deltt / wetfront_cell(t); % / porosity_m(t); 
            (*gv.c_i) =  ( (*gv.c_i) * gv.porosity_i_prev + (*gv.exchange_si) * gv.porosity_s_prev ) / gv.porosity_i; // / porosity_m(t);
            (*gv.c_s) = ( (*gv.c_s) * gv.porosity_s_prev - (*gv.exchange_si) * gv.porosity_s_prev ) / gv.porosity_s;
            //c_s_new_i = c_s_new_i  - exchange ; 

            // Exchange with immobile phase (just exchange)
            (*gv.exchange_im)  = gp.alpha/gv.porosity_m_prev * ((*gv.c_i) - (*gv.c_m)) ; 
            // exchange = alpha .* (c_i_new_i - c_m_new_i) ; 
            // limit the flux to the available material
            
            // limit the flux to the available material
            for(ix=1;ix<gv.nl ;ix++){
                for(iy=1;iy<gv.nh ;iy++){
                    if ((*gv.exchange_im)(ix,iy) > 0){
                        (*gv.exchange_im)(ix,iy) = std::min((*gv.exchange_im)(ix,iy),(*gv.c_i)(ix,iy));
                    }else if((*gv.exchange_im)(ix,iy) < 0){
                        msg = "PROBLEM: negative exchange";
                        print_screen_log(logPULSEfile,msg);
                    }
                }
            };

            (*gv.c_m) =  ( (*gv.c_m) * gv.porosity_m_prev + (*gv.exchange_im) * gv.porosity_i_prev ) / gv.porosity_m; // / porosity_m(t);
            (*gv.c_i) = ( (*gv.c_i) * gv.porosity_i_prev - (*gv.exchange_im) * gv.porosity_i_prev ) / gv.porosity_i;

            // add all immobile and solid slow that melted from the last cell) 
            if(upperboundary_cell_new != upperboundary_cell_prev){
                (*gv.c_m)(arma::span(0,gv.nl),upperboundary_cell_new) = ( (*gv.c_m)(arma::span(0,gv.nl),upperboundary_cell_new) * gv.porosity_m_prev
                    + (*gv.c_m)(arma::span(0,gv.nl),upperboundary_cell_prev) * gv.porosity_m_prev
                    + (*gv.c_s)(arma::span(0,gv.nl),upperboundary_cell_prev) * gv.porosity_s_prev
                    + (*gv.c_i)(arma::span(0,gv.nl),upperboundary_cell_prev) * gv.porosity_i_prev ) / gv.porosity_m; // gv.porosity_m;

                (*gv.c_m)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
                (*gv.c_s)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
                (*gv.c_i)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
            }
            // if porosities are too small, they create instability
        }else{
            gv.porosity_i = 0;
            gv.porosity_m = 1;
            gv.porosity_s = 0;
            (*gv.c_m)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
            (*gv.c_s)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
            (*gv.c_i)(arma::span(0,gv.nl),arma::span(0,upperboundary_cell_new)) = 0;
        }
    
    
        //- removed this code; check if it is important) // c_m_new_i(1,:) = c_m_new_i(3,:);
        //- removed this code; check if it is important) // c_m_new_i(2,:) = c_m_new_i(3,:);   
        //- removed this code; check if it is important) // c_m_new_i(end,:) = c_m_new_i(end-2,:);   
        //- removed this code; check if it is important) // c_m_new_i(end-1,:) = c_m_new_i(end-2,:);  

        //- removed this code; check if it is important) // c_i_new_i(1,:) = c_i_new_i(3,:);
        //- removed this code; check if it is important) // c_i_new_i(2,:) = c_i_new_i(3,:);   
        //- removed this code; check if it is important) // c_i_new_i(end,:) = c_i_new_i(end-2,:);   
        //- removed this code; check if it is important) // c_i_new_i(end-1,:) = c_i_new_i(end-2,:);  

        //- removed this code; check if it is important) // c_s_new_i(1,:) = c_s_new_i(3,:);
        //- removed this code; check if it is important) // c_s_new_i(2,:) = c_s_new_i(3,:);  
        //- removed this code; check if it is important) // c_s_new_i(end,:) = c_s_new_i(end-2,:);   
        //- removed this code; check if it is important) // c_s_new_i(end-1,:) = c_s_new_i(end-2,:);  


        // Save
        //c_m(:,1:wetfront_cell_new) = c_m_new_i;
        //c_i(:,1:wetfront_cell_new) = c_i_new_i;
        //c_s(:,1:wetfront_cell_new) = c_s_new_i;

        // Compute flux
        //c_export= [c_export, c_m(round(nl-2/2),end-1)];
        //cflux = [cflux, c_export(t) * q];



        //- removed this code; check if it is important) //if (t == 2){
        //- removed this code; check if it is important) //    fileID = []; // to save the output file ID for subsequent access to save the following time steps
        //- removed this code; check if it is important) //    tcum_prev = 0;
        //- removed this code; check if it is important) //}

          if (tcum>=print_next)
        {
             end = std::chrono::system_clock::now();
             elapsed_seconds = end-start;
             
              outwritestatus = print_results(gv,gp,std::round(print_next),gp.print_step,elapsed_seconds);
                
            if(outwritestatus == true) 
            {
                std::cout << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(gp.print_step/60) << " || Time elapsed (min): " << elapsed_seconds.count()/60 << std::endl;
                logPULSEfile << "Saved: '" << print_next << ".txt' || time step (min): " << std::to_string(gp.print_step/60) << " || Time elapsed (min): " << std::to_string(elapsed_seconds.count()/60) + "\n";
                print_next += gp.print_step;
                start = std::chrono::system_clock::now();
            } else
            {
                msg = "Problem when saving the results:" + print_next;
                print_screen_log(logPULSEfile,msg);
                abort();
            }
             
         }

    }
    
    
}



int main(int argc, char** argv) 
{   
    
    int H_local,L_local,h_layer,l_layer,nl,nh;
    std::string sim_purp,qmelt_file,msg;
    std::ofstream logPULSEfile ("log.pulse");
    
    // Assign global parameters
    globalpar gp; 
        
    msg = "PULSE: multi-phase multi-layer snowpack chemistry model";
    print_screen_log(logPULSEfile,msg);    
    
    // read simulation setup
    read_simset(gp,sim_purp, &H_local,&L_local,&h_layer,&l_layer,qmelt_file,logPULSEfile);   
    
    // create mesh
    checkmesh(&H_local,&L_local,&h_layer,&l_layer,&nl,&nh,logPULSEfile);
    
    // Asign global variables (heap)
    globalvar gv(nl,nl); 
    (gv.snowH) = H_local;
    (gv.snowL) = L_local;
    (gv.snowh) = h_layer;
    (gv.snowl) = l_layer;
     
    // read snowmelt input
    read_qmelt(gv,qmelt_file,logPULSEfile);
    
    // call the main PULSE model
    PULSEmodel(gp,gv,logPULSEfile);

    // Simulation completed
    msg = "----- Simulation complete -----";
    print_screen_log(logPULSEfile,msg);   
    
    return 0;
}

