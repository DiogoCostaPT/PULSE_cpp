

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
bool print_results(globalvar& gv,globalpar& gp, int print_tag, unsigned int print_step, 
        std::chrono::duration<double> elapsed_seconds,std::string* results_flname)
{

    unsigned int il,ih;
    int a = 0;
    int nh_l = gv.nh;
    
    std::string tprint = *results_flname + "/" + std::to_string(print_tag); 
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
            //filedataR(a,3) = (*gv.c_i).at(il,ih); 
            filedataR(a,3) = (*gv.c_s).at(il,ih); 
            filedataR(a,4) = (*gv.vfrac2d_m).at(il,ih);
            filedataR(a,5) = (*gv.vfrac2d_s).at(il,ih);
            filedataR(a,6) = (*gv.v_liq).at(il,ih)*1000*1000; // m*m*m -> mm*mm*m
            filedataR(a,7) = (*gv.v_swe).at(il,ih)*1000*1000; // m*m*m -> mm*mm*m
            filedataR(a,8) = (*gv.v_air).at(il,ih)*1000*1000; // m*m*m -> mm*mm*m
            //filedataR(a,9) = (*gv.exchange_si).at(il,ih); 
            //filedataR(a,10) = (*gv.exchange_is).at(il,ih); 
            a = a + 1;
        }
    }
   
    arma::mat filedata(std::max(0,a-1),9); 
    if (a>0){
        filedata = filedataR(arma::span(0,std::max(0,a-1)),arma::span(0,8));
    }
    
    int num_export_var = 9;
    arma::field<std::string> header(num_export_var);
    header(0) = "yh [-]";
    header(1) = "xl [-]";
    header(2) = "c_m [user_defined]";
    header(3) = "c_s [user_defined]";
    header(4) = "vfrac_liqwater [-]";
    header(5) = "vfrac_swe [-]";
    header(6) = "v_liq [mm*mm*m]";
    header(7) = "v_swe [mm*mm*m]";
    header(8) = "v_air [mm*mm*m]";     
    
    bool outwritestatus =  filedata.save(arma::csv_name(tprint, header));
    return outwritestatus;
}
