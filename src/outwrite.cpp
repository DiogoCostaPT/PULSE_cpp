

#include "outwrite.h"

/* *****
 * Print to console and log.pulse file 
 ***** */
void print_screen_log(std::ofstream* logPULSEfile,std::string* msg)
{
    //try{
        std::cout << (*msg) << std::endl;
        (*logPULSEfile) << (*msg) + "\n";   
        
    //} catch(const std::exception& e){
    //}
    
}

/* *****
 * Print results   
 * **** */
bool print_results(globalvar& gv,globalpar& gp, int print_tag, unsigned int print_step, std::chrono::duration<double> elapsed_seconds)
{

    unsigned int il,ih;
    int a = 0;
    int nh_l = gv.nh;
    
    std::string tprint = "Results/" + std::to_string(print_tag); 
    std::string filext(".txt");
    tprint += filext;

    arma::mat filedataR(gv.nl*gv.nh,9); 
    
    for(ih=0;ih<gv.nh;ih++)
    {
        for(il=0;il<gv.nl;il++)
        {
            filedataR(a,0) = nh_l - ih - 1;// * gv.snowh;  
            filedataR(a,1) = il;// * gv.snowl;  
            filedataR(a,2) = (*gv.c_m).at(il,ih); 
            filedataR(a,3) = (*gv.c_i).at(il,ih); 
            filedataR(a,4) = (*gv.c_s).at(il,ih); 
            filedataR(a,5) = (gv.vfrac_m); 
            filedataR(a,6) = (gv.vfrac_s); 
            filedataR(a,7) = (*gv.exchange_si).at(il,ih); 
            filedataR(a,8) = (*gv.exchange_im).at(il,ih); 
            a = a + 1;
        }
    }
   
    arma::mat filedata(std::max(0,a-1),8); 
    filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,8));
    
    bool outwritestatus =  filedata.save(tprint,arma::csv_ascii);
    return outwritestatus;
}
