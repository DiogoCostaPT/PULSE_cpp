

folder_loc = '/media/dcosta/data/megasync/ec_main/models/pulse/code/code_matlab_original/Svalbard_Snownet/';


%% Generate 0.txt file 
    
chemistry_file = 'BRG_data.xlsx';
species = 'NO3';

dataraw_chem = xlsread([folder_loc,chemistry_file],species);

time_obs_chem = dataraw_chem(:,1) + 693960;

depths_obs = [10,20,30,40,50,60,70,80,90,100];		

NO3_conc_ppb = dataraw_chem(1,2:end); % ppb
%NO3_conc_ppb = [281.766753953895,112.96520293409,162.480671212122,102.504849889069,116.591016244802,...
%    92.9862893794549,141.758823073791,100.783465395752,92.4981427755126,135.784890366924]; % ppb

NO3_conc_mgl = NO3_conc_ppb / 1000;

snow_L = 2;
snow_h = 0.2; % 5 mm
snow_l = 0.2; % 1 mm

% calc

snow_H = depths_obs(end);

cell_h_num = snow_H/snow_h;
cell_l_num = snow_L/snow_l;

depths_obs_flip = fliplr(depths_obs);

file_0txt = [];

cm_0 = 0;
poros_m = 0.008;
poros_s = 0.991;

for hci = 0:cell_h_num-1
    for lci = 0:cell_l_num-1

    hi = hci * snow_h;

    iloc_max = find(depths_obs_flip>=hi);

    cs_0 = NO3_conc_mgl(iloc_max(end));

    file_0txt = [file_0txt;[hci,lci,cm_0,0,cs_0,poros_m,poros_s,0,0]];

    end

end


%% Generate Qmelt data
   
meteo_file = 'Meteo_2014-2015.xlsx';
dataraw_meteo = importdata([folder_loc,meteo_file]);

time_meteo = dataraw_meteo.data(:,1) + 695422 ;

timestart_chem = time_obs_chem(1);
timeend_chem = time_obs_chem(end);

i_tstart = find(time_meteo==timestart_chem);
i_tend = find(time_meteo==timeend_chem);


dataraw_meteo_relsubset = dataraw_meteo.data(i_tstart:i_tend,:);

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

% PREC_type_ts: 2 is liquid, 3 is snow, 30 and 32 a kind of slush snow
RAIN_ts = PREC_ts;
RAIN_ts.Data(PREC_type_ts.Data~=2) = 0;

SNOWfall_ts = PREC_ts;
SNOWfall_ts.Data(PREC_type_ts.Data==2) = 0; % daily


figure
plot(RAIN_ts)
hold on
plot(SNOWfall_ts)
legend('RAIN_tc','SNOWfall_tc')
grid on
ylabel('mm')


% for model

time_pulse = [0; cumsum(etime(datevec(SNOWfall_ts.Time(2:end)),datevec(SNOWfall_ts.Time(1:end-1))))];
time_pulse = SNOWfall_ts.Time(1:end);
snowaccum_2pulse = SNOWfall_ts.Data/(3600*24); % minus because it is accumulation

Snowfall_file = [time_pulse,snowaccum_2pulse,snowaccum_2pulse*0+0.1];