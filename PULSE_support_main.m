
% General info
masterfile = 'simset.pulse';
pulse_dir = 'bin';
results_dir = 'Results';

% Run pulse once?
Run_pulse_flag = 1;
Clean_results_folder_except_IC_flag = 1;
IC_file = '0.txt';

% Plot results in "Results" directory?
Plot_results_flag = 0;
chemical_species = 'NO3';
col_li = 5; % vertical cell to print results
Obs_file = 'BRG_data.xlsx';

% Sensitivity analysys
Sens_analysys_flag = 1;



%% Run PULSE
if Run_pulse_flag; PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,pulse_dir,...
        results_dir,IC_file,masterfile); end

%% Plot Results
if Plot_results_flag; PULSE_support_plot_results(pulse_dir,results_dir,...
            chemical_species,col_li,masterfile,Obs_file); end








