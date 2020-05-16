

function PULSE_support_Sens_process(Obs_file,chemical_species,col_li,...
       masterfile)

sens_folder = [pwd,'/Sensitivity_analysis/'];
windowtitle = 'Select the folder conaiting the Sensitivity batch runs'; 
sensbatch_dir = uigetdir(sens_folder,windowtitle);

resultsreport_fullpath = [sensbatch_dir,'/Sens_report.txt'];
folder_content = dir(sensbatch_dir);
folder_content = {folder_content.name};
folder_content(~contains(folder_content,'Run_')) = [];
num_runs = numel(folder_content);

% get obs file data
% [depth_fixed_int,depth_corr,Time_data,Obs_data,elev_meas] = PULSE_support_Get_obs_data(Obs_file,...
%     chemical_species);
 
% import Obs_data data
%A_D_all = zeros(num_runs,1)*NaN;
%ALPHA_IE_all = zeros(num_runs,1)*NaN;
%Nash = zeros(num_runs,1)*NaN;
%Bias = zeros(num_runs,1)*NaN;
%RMSE = zeros(num_runs,1)*NaN;

hbar = parfor_progressbar(num_runs, 'Calculating Nash, RMSE and Bias...');
%run_no = str2num(char(run_no));
for run_i = 1:num_runs
    hbar.iterate(1)
    
    sim_i = folder_content{run_i};
    sim_i_folder = [sensbatch_dir,'/',sim_i];
    results_dir = [sim_i_folder,'/Results'];
    resreport_fullpath = [sim_i_folder,'/Sens_results.txt'];
    
    reportexists = exist(resreport_fullpath,'file');
    
    if reportexists ~= 0
        disp(['SensProcess: simulation already processed (',sim_i,')'])
        continue;
    end
    
    % get parameters
    param_folder = [sim_i_folder,'/Sens_info.txt'];
    param = importdata(param_folder);
    A_D = param.data(1);
    ALPHA_IE = param.data(2);
    
    % load model results
    [time_sim_elapsec,h_layers_max,c_m,c_s,c_total,poros_m,poros_s] = ...
         PULSE_support_load_pulse_results(results_dir,col_li);
     
    % Retrieve Obs data for c_total
    [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec,...
        cmax_ctotal] = PULSE_support_GetTrans_obs_data(Obs_file,c_total,...
        chemical_species); 
     
    % get comment and timenum from masterfile
    [comment,time_sim,H_LAY,~] = PULSE_support_Getinfo_masterfile(...
            time_sim_elapsec,masterfile);
    H_LAY = H_LAY/10; % mm to cm
    ih = 0:H_LAY:(h_layers_max-1)*H_LAY;
    
     % generate mesh grid for plot surf
    [Hmesh,Tmesh] = meshgrid(ih,time_sim);
             
    %Model_data_interc = interp1(Data_time(1:end-1),Data(1:end-1),T_WQsort);
    Model_data_interc = interp2(Hmesh,Tmesh,c_total,Y_obs_mesh,X_obs_mesh);
    
    % remove NaNs
    nanloc = find(isnan(Model_data_interc));
    Model_data_interc(nanloc) = [];
    Z_obs_mesh(nanloc) = [];
    
    nanloc = find(isnan(Z_obs_mesh));
    Model_data_interc(nanloc) = [];
    Z_obs_mesh(nanloc) = [];
    
    %figure
    %scatter(Z_obs_mesh,Model_data_interc)
    %limmin = min(min(Model_data_interc));
    %limmax = max(max(Model_data_interc));
    %hold on
    %plot([limmin limmin],[limmax limmax],'k')
    %xlim([limmin limmax])
    %grid on
    %xlabel('Obs (mg/l)')
    %ylabel('Model (mg/l)')
    
    Model_data_interc_use = Model_data_interc;
    WQ_use = Z_obs_mesh;
    
    % Nash
    numerator=(WQ_use-Model_data_interc_use).^2;
    denominator=(WQ_use-mean(WQ_use)).^2;
    Nash =1-(sum(numerator)/sum(denominator));
    
    % RMSE
    Sumcal = (Model_data_interc_use-WQ_use).^2;
    numerator = sum(Sumcal);
    n=numel(WQ_use);
    RMSE=(numerator/n)^(1/2);
    
    % BIAS
    numerator = sum(WQ_use);
    denominator = sum(Model_data_interc_use);
    Bias = numerator/denominator-1;
    
    printcell = {'A_D,','ALPHA_IE','Nash','Bias','RMSE';
                  A_D,ALPHA_IE,Nash,Bias,RMSE};
    writecell(printcell,resreport_fullpath,'delimiter',',')
      
    
end
 hbar.close();

% import Obs_data data
ResSens_all = cell(num_runs,5);
ResSens_all(1,1:end) = {'A_D,','ALPHA_IE','Nash','Bias','RMSE'};

% Save now in one file
hbar = parfor_progressbar(num_runs, 'Calculating Nash, RMSE and Bias...');
%run_no = str2num(char(run_no));
parfor run_i = 1:num_runs
    hbar.iterate(1)
    
    sim_i = folder_content{run_i};
    sim_i_folder = [sensbatch_dir,'/',sim_i];
    resreport_fullpath = [sim_i_folder,'/Sens_results.txt'];
    
    reportexists = exist(resreport_fullpath,'file');
    
     if reportexists ~= 0
        disp(['SensProcess: cannot find sens report for (',sim_i,') -> skipped'])
        continue;
     else
        dataall = importdata(resreport_fullpath);
        ResSens_all(run_i,:) = dataall.data;
     end
 
end
 hbar.close();
 
 writecell(ResSens_all,resultsreport_fullpath,'delimiter',',');