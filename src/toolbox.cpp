

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
   std::vector<char*> filenames; //stringvec filenames, filename_i;
   const char *filename_i;
   char *simnum_str_i;
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
        simnum_str_i = (char*) malloc(sizeof(filename_i)-2);
        strncpy (simnum_str_i, filename_i, sizeof(filename_i)-2);
        simnum = atoi(simnum_str_i);
        timestart = std::max(timestart,simnum);
        free(simnum_str_i);
   }
   
   free(entry);
   return timestart;
}

/* 
 * Identify the snowpack mesh based on file 0.txt 
 */
void checkmesh2(int* H_local,int* L_local,int* h_layer,int* l_layer,int* nh,int* nl,std::ofstream* logPULSEfile)
{
    unsigned int a, timstart;
    
    arma::mat filedata; 
    std::string init_file, msg;
    
    timstart = findLastStep("Results/"); // list the results files to get the last time step
    init_file = "Results/" + std::to_string(int(timstart)) + ".txt";   
    bool flstatus = filedata.load(init_file,arma::csv_ascii);
    
    *nh = 0;
    *nl = 0;

    if(flstatus == true) 
    {
        for(a=0;a<filedata.col(1).n_elem;a++)
        {
            (*nh) = std::max((*nh), int(filedata(a,0)) + 1);  
            (*nl) = std::max((*nl), int(filedata(a,1)) + 1);  
          
        } 
        (*H_local) =  (*nh) * (*h_layer);
        (*L_local) =  (*nl) * (*l_layer);
        
        msg = "Mesh: identified ";
        print_screen_log(logPULSEfile,&msg);  
        
    }else{
        msg = "Mesh: not identified -> Results/*.txt file(s) missing";
        print_screen_log(logPULSEfile,&msg);  
        std::abort();
    }     
}

//void checkmesh(int* H_local,int* L_local,int* h_layer,int* l_layer,int* nh,int* nl,std::ofstream* logPULSEfile)
//{
//
//    std::string msg;
//    
//    double a = double(*H_local)/double(*h_layer);
//    double b = double(*L_local)/double(*l_layer);
//    
//    if(floor(a)==ceil(a) && floor(b)==ceil(b)){      
//        (*nh) = a;
//        (*nl) = b;
//        msg = "Snowpack mesh: created";
//        print_screen_log(logPULSEfile,&msg);
//    }else{
//        if (floor(a)==ceil(a)){
//        msg = "Snowpack mesh: H = " + std::to_string ((*H_local)) + 
//                " mm (snowpack depth) is not divisible by h_layer = " + std::to_string((*h_layer)) +
//                " mm (grid thickness) -> change the simulation setting in file 'simset.pulse'";
//        print_screen_log(logPULSEfile,&msg);
//        }
//        if (floor(b)==ceil(b)){
//        msg = "Snowpack mesh: L = " + std::to_string ((*L_local)) + 
//                " mm (snowpack horizontal length) is not divisible by l_layer = " + std::to_string((*l_layer)) +
//                "mm (grid horizontal length) -> change the simulation setting in file 'simset.pulse'";
//        print_screen_log(logPULSEfile,&msg);
//        }
//        std::abort();
//    }
//    
//}