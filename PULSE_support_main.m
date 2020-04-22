
% General info
masterfile = 'simset';
pulse_dir = 'bin';
results_dir = 'bin/Results';
Obs_file = 'BRG_data.xlsx';
chemical_species = 'NO3';
col_li = 5; % vertical cell to print results

% Run pulse once?
Run_pulse_flag = 0;
Clean_results_folder_except_IC_flag = 1;
IC_file = '0.txt';

% Plot results in "Results" directory?
Plot_results_flag = 0;

% Sensitivity analysys run?
Sens_run_flag = 1;
num_samples = 200;
A_D_max = 0.0003;           
ALPHA_IE_max = 0.000003;

% Analyse the runs in Sensitivity_analysis folder
Sens_analysis_flag = 1;

% Plot Sens results
Sens_plot_results_flag = 1;



%% Run PULSE (once)
if Run_pulse_flag; PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,pulse_dir,...
        results_dir,IC_file,masterfile); end

%% Plot Results
if Plot_results_flag; PULSE_support_plot_results(pulse_dir,results_dir,...
            chemical_species,col_li,masterfile,Obs_file); end
    
%% Sensitivity runs (multiple runs)
if Sens_run_flag; PULSE_support_Sens_analysys_run(masterfile,pulse_dir,num_samples,...
        A_D_max,ALPHA_IE_max); end

% Sensitivity analysis (processing the results in Sensitivity_analysis
if Sens_analysis_flag; PULSE_support_Sens_process(Obs_file,chemical_species,col_li,...
        masterfile); end

% Sensitivity results plotting
if Sens_plot_results_flag ;PULSE_support_Sens_plot() ; end








