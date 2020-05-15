
% general info
folder_loc = '/media/dcosta/data/megasync/ec_main/models/pulse/code/code_matlab_original/Svalbard_Snownet/';
meteo_file = 'Meteo_2014-2015.xlsx';
chemistry_file = 'BRG_data.xlsx';
species = 'NO3';

% IC file
gen_0txtfile_flag = 0;   
snow_L = 100; % mm
snow_H = 1000; % 100 cm * 10 = 1000 mm
snow_h = 10; % 10 mm
snow_l = 10; % 10 mm
v_frac_air_init = 0.05; % volume percertage

% qmelt and meteo files
gen_prec_and_qmelt_T_index_files_flag = 1;
snowmelt_method = 1; % 0)T-index, 1) CRHM output
T_index_coef = 10; % only used if snowmelt_method = 0
crhmoutput_dir = '/media/dcosta/data/megasync/ec_main/models/crhm/support/PROJECTS/Svalbard/CRHM_output_1.txt'; % only used if snowmelt_method =  1;
new_qmeltfile_name = 'qcmelt_svalbard_test_crhm';
new_meteofile_name = 'meteo_svalbard_test';


% Densities (as in PULSE)
rho_s = 917;  % kg.m-3 at 0 degrees
rho_m = 998.8; % kg.m-3 at 0 degrees
rho_frshsnow_init = 320;


%% Generate 0.txt file 
if gen_0txtfile_flag == 1
  
    dataraw_chem = xlsread([folder_loc,chemistry_file],species);

    time_obs_chem = dataraw_chem(:,1) + 693960;

    depths_obs = [10,20,30,40,50,60,70,80,90,100] * 10; % cm -> mm

    NO3_conc_ppb = dataraw_chem(1,2:end); % ppb
    %NO3_conc_ppb = [281.766753953895,112.96520293409,162.480671212122,102.504849889069,116.591016244802,...
    %    92.9862893794549,141.758823073791,100.783465395752,92.4981427755126,135.784890366924]; % ppb

    NO3_conc_mgl = NO3_conc_ppb / 1000; %mg/l

    % calc

    %snow_H = depths_obs(end);

    cell_h_num = snow_H/snow_h;
    cell_l_num = snow_L/snow_l;

    depths_obs_flip = fliplr(depths_obs);

    file_0txt = [];

    cm_0 = 0;
    v_liqwater = 0; % m3
    v_swe = snow_h * snow_l * rho_frshsnow_init/rho_m; % m3
    v_air = snow_h * snow_l * v_frac_air_init;
    
    vfrac_s = v_swe / (v_swe + v_liqwater);
    vfrac_m = v_liqwater / (v_swe + v_liqwater);

    for hci = 0:cell_h_num-1
        for lci = 0:cell_l_num-1

        hi = hci * snow_h;

        iloc_max = find(depths_obs_flip>=hi);

        cs_0 = NO3_conc_mgl(iloc_max(end));

        file_0txt = [file_0txt;[hci,lci,cm_0,cs_0,vfrac_m,vfrac_s,v_liqwater,v_swe,v_air]];

        end

    end
    
    file_0txt_cell = num2cell(file_0txt);
    header = {};
    header{1} = 'yh [-]';
    header{2} = "xl [-]";
    header{3} = "c_m [user_defined]";
    header{4} = "c_s [user_defined]";
    header{5} = "vfrac_liqwater [-]";
    header{6} = "vfrac_swe [-]";
    header{7} = "v_liqwater [mm*mm*m]";
    header{8} = "v_swe [mm*mm*m]";
    header{9} = "v_air [mm*mm*m]";     
    
    file_0txt_cell = [header;file_0txt_cell];
    
    writecell(file_0txt_cell,'0.txt','Delimiter',',')

end
%% Generate Qmelt data
if gen_prec_and_qmelt_T_index_files_flag == 1
    
    % meteo data
    dataraw_meteo = importdata([folder_loc,meteo_file]);
    time_meteo = dataraw_meteo.data(:,1) + 695422 ;
    
    % chem data
    dataraw_chem = xlsread([folder_loc,chemistry_file],species);
    depths_obs = [10,20,30,40,50,60,70,80,90,100];
    time_obs_chem = dataraw_chem(:,1) + 693960;
    timestart_chem = time_obs_chem(1);
    timeend_chem = time_obs_chem(end);
    
    % sub-set meteo data to chem data time period
    i_tstart = find(time_meteo==timestart_chem);
    i_tend = find(time_meteo==timeend_chem);
    dataraw_meteo_relsubset = dataraw_meteo.data(i_tstart:i_tend,:);

    % prepare data
    TIME = dataraw_meteo_relsubset(:,1) + 695422;
    TEMP = dataraw_meteo_relsubset(:,3);
    RH = dataraw_meteo_relsubset(:,4);
    PRESS = dataraw_meteo_relsubset(:,5);
    WS = dataraw_meteo_relsubset(:,6);
    WD = dataraw_meteo_relsubset(:,7);
    RAD = dataraw_meteo_relsubset(:,8);
    PREC = dataraw_meteo_relsubset(:,9);
    PREC_type = dataraw_meteo_relsubset(:,10);

    TEMP_ts = timeseries(TEMP,TIME);
    RH_ts = timeseries(RH,TIME);
    PRESS_ts = timeseries(PRESS,TIME);
    WS_ts = timeseries(WS,TIME);
    WD_ts = timeseries(WD,TIME);
    RAD_ts = timeseries(RAD,TIME);
    PREC_ts = timeseries(PREC,TIME);
    PREC_type_ts = timeseries(PREC_type,TIME);

    % plot all data
    figure
    subplot(3,3,1)
    plot(TEMP_ts)
    title('TEMP_ts')
    grid on
    subplot(3,3,2)
    plot(RH_ts)
    title('RH_ts')
    grid on
    subplot(3,3,3)
    plot(PRESS_ts)
    title('PRESS_ts')
    grid on
    subplot(3,3,4)
    plot(WS_ts)
    title('WS_ts')
    grid on
    subplot(3,3,5)
    plot(WD_ts)
    title('WD_ts')
    grid on
    subplot(3,3,6)
    plot(RAD_ts)
    title('RAD_ts')
    grid on
    subplot(3,3,7)
    plot(PREC_ts)
    title('PREC_ts')
    grid on
    subplot(3,3,8)
    plot(PREC_type_ts)
    title('PREC_type_ts')
    grid on

    %% PREC_type_ts: 2 is liquid, 3 is snow, 30 and 32 a kind of slush snow
    RAIN_ts = PREC_ts;
    RAIN_ts.Data(PREC_type_ts.Data~=2) = 0;

    SNOWfall_ts = PREC_ts;
    SNOWfall_ts.Data(PREC_type_ts.Data==2) = 0; % daily
    
           
    % plot
    figure
    subplot(1,2,1)
    plot(RAIN_ts)
    hold on
    plot(SNOWfall_ts)
    legend('RAIN_tc','SNOWfall_tc')
    grid on
    ylabel('mm')
    subplot(1,2,2)
    plot(cumsum(RAIN_ts.data))
    hold on
    plot(cumsum(SNOWfall_ts.data))
    legend('RAIN_tc','SNOWfall_tc')
    grid on
    ylabel('mm')

    % data to copy past to model's meteo file (time and prec vol)
    time_pulse = SNOWfall_ts.Time(1:end);
    snowaccum_2pulse = SNOWfall_ts.Data; % mm acumul
    rainaccum_2pulse = RAIN_ts.Data;
    temp_2pulse = TEMP_ts.Data;

    % now chem data
    data_toplayer = dataraw_chem(:,2)/1000; % ppb -> ppm
    snowprec_chem = interp1(time_obs_chem,data_toplayer,time_pulse);
    
    % join
    time_pulse_sec = [0; cumsum(etime(datevec(SNOWfall_ts.Time(2:end)),datevec(SNOWfall_ts.Time(1:end-1))))];
    snowaccum_2pulse(isnan(snowaccum_2pulse)) = 0;
    snowprec_chem(isnan(snowprec_chem)) = 0;
    Meteo_file = [time_pulse_sec,...
                  temp_2pulse,...
                  rainaccum_2pulse,...
                  snowaccum_2pulse,...
                  snowprec_chem];
    

    Meteo_file_cell = num2cell(Meteo_file);
    header = {};
    header{1} = 'time (sec)';
    header{2} = "temperature [degree celsius]";
    header{3} = "rain [mm/deltatime]";
    header{4} = "snowfall [mm/deltatime]";
    header{5} = "prec_conc [user_defined]";

    Meteo_file_cell = [header;Meteo_file_cell];
    
    writecell(Meteo_file_cell,new_meteofile_name,'Delimiter',',')

        
    %% QMELT
    if snowmelt_method == 0
        qmelt_estim_ts = max(TEMP_ts.data * T_index_coef,0);
        qmelt_estim_ts(isnan(qmelt_estim_ts)) = 0;
        qmelt_file = [time_pulse_sec,qmelt_estim_ts];
    else
       crhm_output = readtable(crhmoutput_dir);     
       timecrhm = str2double(crhm_output.time(2:end)) + 693960;
       snowmelt_data = str2double(crhm_output.snowmelt_int_1_);
       
       date_start = datenum('25-03-2015','dd-mm-yyyy');
       date_end = datenum(date_start + seconds(5616000));
       
       date_start_loc = find(timecrhm == date_start);
       date_end_loc = find(timecrhm == date_end);
       
       time_crhm_relvt = timecrhm(date_start_loc:date_end_loc);
       time_crhm_relvt_secdiff = etime(datevec(time_crhm_relvt(2:end)),...
                                datevec(time_crhm_relvt(1:end-1)));
       time_crhm_relvt_secdiff_cumsum = [0;cumsum(time_crhm_relvt_secdiff)];
       snowmelt_data_relvt = snowmelt_data(date_start_loc:date_end_loc);
       
       qmelt_file = [time_crhm_relvt_secdiff_cumsum,snowmelt_data_relvt];
       
    end
    
    qmelt_file_cell = num2cell(qmelt_file);
    header = {};
    header{1} = 'time (sec)';
    header{2} = "qmelt [mm/deltatime]";

    qmelt_file_cell = [header;qmelt_file_cell];
    
    writecell(qmelt_file_cell,new_qmeltfile_name,'Delimiter',',')
    
     % plot
    figure
    subplot(1,2,1)
    plot(SNOWfall_ts)
    hold on
    qmelt_file_ts = timeseries(qmelt_file(:,2),qmelt_file(:,1));
    plot(qmelt_file_ts)
    legend('SNOWfall_tc','Qmelt')
    grid on
    ylabel('mm')
    subplot(1,2,2)
     plot(cumsum(SNOWfall_ts.data))
    hold on
    plot(cumsum(qmelt_file_ts.data))
    legend('SNOWfall_tc','Qmelt')
    grid on
    ylabel('mm')

end

