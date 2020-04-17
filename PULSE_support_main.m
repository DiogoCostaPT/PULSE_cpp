
% General info
masterfile = 'simset.pulse';
pulse_dir = 'bin';
results_dir = 'bin/Results';

% Run pulse once?
Run_pulse_flag = 1;
Clean_results_folder_except_IC_flag = 1;
IC_file = '0.txt';

% Plot results in "Results" directory?
Plot_results_flag = 0;
chemical_species = 'NO3';
col_li = 5; % vertical cell to print results
Obs_file = 'BRG_data.xlsx';

% Sensitivity analysys?
Sens_analysys__run_flag = 0;
num_samples = 200;
A_D_max = 0.3;           
ALPHA_IE_max = 0.3;


%% Run PULSE (once)
if Run_pulse_flag; PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,pulse_dir,...
        results_dir,IC_file,masterfile); end

%% Plot Results
if Plot_results_flag; PULSE_support_plot_results(pulse_dir,results_dir,...
            chemical_species,col_li,masterfile,Obs_file); end
    
%% Sensitivity analysis (multiple runs)
if Sens_analysys__run_flag; PULSE_support_Sens_analysys_run(masterfile,pulse_dir,num_samples,...
        A_D_max,ALPHA_IE_max); end








