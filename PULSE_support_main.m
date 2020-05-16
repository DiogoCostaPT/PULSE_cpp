
%% GUI
varargout = PULSE_GUI_App();
while(varargout.RUNButton.Value ~=1) 
    pause(0.2);
end

%% Set user data
%masterfile_fullpath = varargout.MasterfilenameEditField.Value; % modset
%pulse_dir = varargout.DirectoryofPULSE_cppEditField.Value; % bin/
%folder_2save_ICfile = varargout.ResultsdirectoryEditField.Value; %'bin/Results';
%chemistry_file = varargout.ObservationfilefullpathEditField.Value; %'BRG_data.xlsx';
%chemical_species = varargout.ChemicalspeciesEditField.Value; %'NO3';
%col_li = 5; % vertical cell to print results

%col_li = round(L_SNOWPACK/L_LAY/2);

% Run pulse once or Sensitivity analysis
optionname = varargout.Switch.Value;
Run_pulse_flag = 0;
Sens_run_flag = 0;
Clean_results_folder_except_IC_flag =0;
Plot_results_flag = 0;
Sens_analysis_flag = 0;
Sens_plot_results_flag = 0;

if strcmp(optionname,'Single RUN')
    Run_pulse_flag = varargout.runsimulationCheckBox.Value;
    Clean_results_folder_except_IC_flag = varargout.deleteResultsfolderexceptICfileCheckBox.Value;
    Plot_results_flag = varargout.PlotresultssinglerunCheckBox.Value;
else
   Sens_run_flag = varargout.RuntestsCheckBox.Value;
   Sens_analysis_flag = varargout.AnalyzeresultsCheckBox.Value; %0;
   Sens_plot_results_flag = varargout.PlotSensresultsCheckBox.Value;
end

%Run_pulse_flag = 1;
IC_file = varargout.edit3.Value; % 0.txt

%Sens_run_flag = 0;
num_samples = varargout.scenariosEditField.Value; %500;
A_D_max = varargout.A_D_maxEditField.Value; %0.0003;           
ALPHA_IE_max = varargout.ALPHA_IE_maxEditField.Value; %0.000003;

% close GUI
close(varargout.figure1)


%% Run PULSE (once)
if Run_pulse_flag; PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,pulse_dir,...
        folder_2save_ICfile,IC_file,masterfile_fullpath); end

%% Plot Results
if Plot_results_flag; PULSE_support_plot_results(pulse_dir,folder_2save_ICfile,...
            chemical_species,col_li,masterfile_fullpath,chemistry_file); end
    
%% Sensitivity runs (multiple runs)
if Sens_run_flag; PULSE_support_Sens_analysys_run(masterfile_fullpath,pulse_dir,num_samples,...
        A_D_max,ALPHA_IE_max); end

% Sensitivity analysis (processing the results in Sensitivity_analysis
if Sens_analysis_flag; PULSE_support_Sens_process(chemistry_file,chemical_species,col_li,...
        masterfile_fullpath); end

% Sensitivity results plotting
if Sens_plot_results_flag ;PULSE_support_Sens_plot() ; end








