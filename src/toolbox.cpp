
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

#include "toolbox.h"

std::string SplitFilename (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(0,found);

}

/* *****
 * Read file names in Results directory 
 ***** */
int findLastStep(const char *path) {

   struct dirent *entry;
   int i, timestart, filenum = 0, simnum;
   std::vector<std::string> filenames; //stringvec filenames, filename_i;
   std::string filename_i;
   //char simnum_str_i;
   DIR *dir = opendir(path);
   
   if (dir != NULL) {
        while ((entry = readdir(dir)) != NULL) {
        filenames.push_back(entry->d_name); // storing the file names
        filenum = filenum + 1;
        }
   }
   closedir(dir);
   
   timestart = 0;
   for(i=2;i<filenum;i++){
       filename_i = filenames[i]; //.assign(filenames[i]); //strcpy(filename_i,(char *).at(&filenames[i]));
        //simnum_str_i = (char*) malloc(sizeof(filename_i)-2);
        //strncpy (&simnum_str_i, &filename_i, sizeof(&filename_i)-2);
        try{
            simnum = std::stoi(filename_i.substr(0,sizeof(filename_i)-4));
            timestart = std::max(timestart,simnum);
        } catch(const std::exception& e){
        }
        //free(simnum_str_i);
   }
   
   free(entry);
   return timestart;
}

/* 
 * Identify the snowpack mesh based on file 0.txt 
 */
void checkmesh2(double* H_local,double* L_local,double* h_layer,double* l_layer,
                int* nh,int* nl,std::ofstream* logPULSEfile,std::string* results_flname)
{
    unsigned int a, timstart;
    
    arma::mat filedata; 
    std::string init_file, msg;
    
    char results_flname_char[(*results_flname).size()+1];
    strcpy(results_flname_char,(*results_flname).c_str());

    timstart = findLastStep(results_flname_char); // list the results files to get the last time step
    init_file = *results_flname + '/' + std::to_string(int(timstart)) + ".txt";   
    bool flstatus = filedata.load(init_file,arma::csv_ascii);
    
    *nh = 0;
    *nl = 0;

    msg = "> Checking MESH... ";
    print_screen_log(logPULSEfile,&msg);  

    if(flstatus == true) 
    {
        for(a=0;a<filedata.n_rows;a++)
        {
            (*nh) = std::max((*nh), int(filedata(a,0)) + 1);  
            (*nl) = std::max((*nl), int(filedata(a,1)) + 1);  
          
        } 
        (*H_local) =  (*nh) * (*h_layer);
        (*L_local) =  (*nl) * (*l_layer);
        
        msg = "    > MESH identified";
        print_screen_log(logPULSEfile,&msg);  
        
    }else{
        msg = "    > MESH fail: Results/*.txt file(s) missing";
        print_screen_log(logPULSEfile,&msg);  
        std::abort();
    }     
}

// Read general matrix file for SNOWPACK = external
bool read_matrixes_ext(globalpar& gp,globalvar& gv,
            std::string* time_file, std::string* v_ice_file,std::string* v_liquid_file,
            std::string* v_ice2liq_1_file,std::string* v_ice2liq_2_file, std::string*fluxQ_file, 
            std::string* prec_c_ext_file,std::ofstream* logPULSEfile)
{
    bool err_flag = false;
    bool flstatus;
    std::string* file_failed, msg;

    msg = "> Preparing simulation... ";
    print_screen_log(logPULSEfile,&msg);  

    if (err_flag == false) flstatus = (*gv.time_ext).load(*time_file,arma::csv_ascii);
    if (flstatus == false) file_failed = time_file;

    if (err_flag == false) flstatus = (*gv.preci_c_ext).load(*prec_c_ext_file,arma::csv_ascii);
    if (flstatus == false) file_failed = prec_c_ext_file;
    
    if (err_flag == false) flstatus = (*gv.v_swe_ext).load(*v_ice_file,arma::csv_ascii);
    if (flstatus == false) file_failed = v_ice_file;
    
    if (err_flag == false) flstatus = (*gv.v_liq_ext).load(*v_liquid_file,arma::csv_ascii);
    if (flstatus == false) file_failed = v_liquid_file;

    if (err_flag == false) flstatus = (*gv.v_ice2liq_1_ext).load(*v_ice2liq_1_file,arma::csv_ascii);
    if (flstatus == false) file_failed = v_ice2liq_1_file;

    if (err_flag == false) flstatus = (*gv.v_ice2liq_2_ext).load(*v_ice2liq_2_file,arma::csv_ascii);
    if (flstatus == false) file_failed = v_ice2liq_2_file;

    if (err_flag == false) flstatus = (*gv.fluxQ_ext).load(*fluxQ_file,arma::csv_ascii);
    if (flstatus == false) file_failed = fluxQ_file;

    if (err_flag == false)
    {
        std::string msg = "    > Completed successfully";
        print_screen_log(logPULSEfile,&msg); 
    }else{
        std::string msg = "    > ERROR reading file: " + (*file_failed);
        print_screen_log(logPULSEfile,&msg); 
        err_flag = true;
        return err_flag;
    }

}

// Function to remove all spaces from a given string 
std::string removeSpaces(std::string str)  
{ 
    str.erase(remove(str.begin(), str.end(), ' '), str.end()); 
    return str; 
} 