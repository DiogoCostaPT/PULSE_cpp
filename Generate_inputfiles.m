

varargout = Generate_inputfiles_GUI();
winopen = 1;
pause(2);

while(winopen == 1)
    
    genetate_masterfile = varargout.GENERATEMASTERFILEButton_2.Value;
    gen_0txtfile_flag = varargout.GENERATEICFILEButton_2.Value;   
    gen_meteo_file_flag = varargout.GENERATEMETEOFILEButton_3.Value;
    gen_qmelt_file_flag = varargout.GENERATEQMELTFILEButton.Value;
    
    if ~genetate_masterfile && ~gen_0txtfile_flag && ~gen_meteo_file_flag && ~gen_qmelt_file_flag  
        pause(2)
        continue; 
    end
    
    % Master file
    masterfile_fullpath = varargout.MASTERFILENAMEfullorrelativepathEditField.Value;
    COMMENT = varargout.HYDRO_SOLVEREditField_2.Value;
    START_TIME = varargout.START_TIMEEditField.Value;
    END_TIME = varargout.END_TIMEEditField.Value;
    PRINT_STEP = varargout.PRINT_STEPsecEditField.Value;
    L_LAY = varargout.L_LAYmmEditField_2.Value; % mm
    H_LAY = varargout.H_LAYmmEditField_2.Value; % 100 cm * 10 = 1000 mm
    VFRAC_AIR_FRESHSNOW = varargout.VFRAC_AIR_FRESHSNOWEditField.Value;
    DENSITY_ICE = varargout.DENSITY_ICEkgm3EditField.Value;
    DENSITY_WATER = varargout.DENSITY_WATERkgm3EditField.Value;
    DENSITY_FRESHSNOW = varargout.DENSITY_FRESHSNOWkgm3EditField.Value;
    A_D = varargout.A_Dm2sEditField.Value;
    ALPHA_IE = varargout.ALPHA_IEEditField.Value;
    COMPFACTOR = varargout.COMPFACTOREditField.Value;
    HYDRO_SOLVER = varargout.HYDRO_SOLVEREditField.Value;
    METEO_FILE = varargout.METEO_FILEfullorrelativepathEditField.Value;
    QMELT_FILE = varargout.QMELT_FILEfullorrelativepathEditField.Value;
    
    % IC file 
    folder_2save_ICfile = varargout.SAVEFILEINFOLDERfullorrelativepathEditField.Value;
    meteo_file = varargout.MeteoObservationsfileEditField.Value;
    chemistry_file = varargout.SnowChemistryfilefullorrelativepathEditField.Value;
    species = varargout.ChemicalspeciesworksheetnameinthefileaboveEditField.Value; 
    H_SNOWPACK = varargout.H_SNOWPACKmmEditField.Value; % 10 mm
    L_SNOWPACK = varargout.L_SNOWPACKmmEditField.Value; % 10 mm
    v_frac_air_init = varargout.VFRAC_AIR_FRESHSNOWEditField.Value; % volume percertage

    % qmelt and meteo files
    snowmelt_method = varargout.snowmeltcalcmethod1Tindex2CRHMEditField_3.Value; % 0)T-index, 1) CRHM output
    T_index_coef = varargout.T_index_coefifsnowmeltcalcmethod1EditField_3.Value; % only used if snowmelt_method = 0
    crhmoutput_dir = varargout.CRHMoutputfilefullorrelativepathEditField.Value; % only used if snowmelt_method =  1;
   
    % reset uneditable boxes
    varargout.SnowChemistryfilefullorrelativepathEditField_2.Value = varargout.SnowChemistryfilefullorrelativepathEditField.Value;
    varargout.ChemicalspeciesworksheetnameinthefileaboveEditField_2.Value = varargout.ChemicalspeciesworksheetnameinthefileaboveEditField.Value; 
    varargout.METEO_FILEfullorrelativepathEditField_2.Value = varargout.METEO_FILEfullorrelativepathEditField.Value;
    varargout.QMELT_FILEfullorrelativepathEditField_2.Value = varargout.QMELT_FILEfullorrelativepathEditField.Value;
    varargout.H_LAYmmEditField_2.Value = varargout.H_LAYmmEditField.Value
    varargout.L_LAYmmEditField_2.Value = varargout.L_LAYmmEditField.Value
    varargout.VFRAC_AIR_FRESHSNOWEditField_2.Value = varargout.VFRAC_AIR_FRESHSNOWEditField.Value;
    
    if genetate_masterfile == 1
        
        varargout.newfilehasbeengeneratedLamp.Color = 'white';
        pause(0.1)
        
        masterfile_txt = {'COMMNET', COMMENT;...
                      'START_TIME', START_TIME;...
                      'END_TIME', END_TIME;...
                      'PRINT_STEP', PRINT_STEP;...
                      'H_LAY_mm', H_LAY;...
                      'L_LAY_mm', L_LAY;...
                      'VFRAC_AIR_FRESHSNOW', VFRAC_AIR_FRESHSNOW;...
                      'DENSITY_ICE', DENSITY_ICE;...
                      'DENSITY_WATER', DENSITY_WATER;...
                      'DENSITY_FRESHSNOW', DENSITY_FRESHSNOW;...
                      'A_D', A_D;...
                      'ALPHA_IE', ALPHA_IE;...
                      'COMPFACTOR', COMPFACTOR;...
                      'QMELT_FILE', QMELT_FILE;...
                      'METEO_FILE', METEO_FILE;...
                      'HYDRO_SOLVER', HYDRO_SOLVER};
                  
       writecell(masterfile_txt,masterfile_fullpath,'Delimiter',' ')
        
        varargout.GENERATEMASTERFILEButton_2.Value = 0;
        varargout.newfilehasbeengeneratedLamp.Color = 'green';
        
    end
    
    %% Generate 0.txt file 
    if gen_0txtfile_flag == 1
        
        varargout.newfilehasbeengeneratedLamp_2.Color = 'white';
        pause(0.1)

        dataraw_chem = xlsread(chemistry_file,species);

        time_obs_chem = dataraw_chem(:,1) + 693960;

        depths_obs = [10,20,30,40,50,60,70,80,90,100] * 10; % cm -> mm

        NO3_conc_ppb = dataraw_chem(2,2:end); % ppb
        %NO3_conc_ppb = [281.766753953895,112.96520293409,162.480671212122,102.504849889069,116.591016244802,...
        %    92.9862893794549,141.758823073791,100.783465395752,92.4981427755126,135.784890366924]; % ppb

        NO3_conc_mgl = NO3_conc_ppb / 1000; %mg/l

        % calc

        %snow_H = depths_obs(end);

        cell_h_num = H_SNOWPACK/H_LAY;
        cell_l_num = L_SNOWPACK/L_LAY;

        depths_obs_flip = fliplr(depths_obs);

        file_0txt = [];

        cm_0 = 0;
        v_liqwater = 0; % m3
        v_swe = H_LAY * L_LAY * DENSITY_FRESHSNOW/DENSITY_WATER; % m3
        v_air = H_LAY * L_LAY * v_frac_air_init;

        vfrac_s = v_swe / (v_swe + v_liqwater);
        vfrac_m = v_liqwater / (v_swe + v_liqwater);

        for hci = 0:cell_h_num-1
            for lci = 0:cell_l_num-1

            hi = hci * H_LAY;

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
        
        folder0txt_path = [folder_2save_ICfile,'/0.txt'];
        folder0txt_path = strrep(folder0txt_path,'//','/');

        writecell(file_0txt_cell,folder0txt_path,'Delimiter',',')
        
        varargout.GENERATEICFILEButton_2.Value = 0;
        varargout.newfilehasbeengeneratedLamp_2.Color = 'green';
        

    end
    %% Generate Qmelt data
    if gen_meteo_file_flag == 1 || gen_qmelt_file_flag == 1

        % meteo data
        dataraw_meteo = importdata(meteo_file);
        time_meteo = dataraw_meteo.data(:,1) + 695422 ;

        % chem data
        dataraw_chem = xlsread(chemistry_file,species);
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

       %{
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
        %}

        %% PREC_type_ts: 2 is liquid, 3 is snow, 30 and 32 a kind of slush snow
        RAIN_ts = PREC_ts;
        RAIN_ts.Data(PREC_type_ts.Data~=2) = 0;

        SNOWfall_ts = PREC_ts;
        SNOWfall_ts.Data(PREC_type_ts.Data==2) = 0; % daily

        %{
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
        %}

        if gen_meteo_file_flag == 1
            
            varargout.newfilehasbeengeneratedLamp_3.Color = 'white';
            pause(0.1)
            
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

            writecell(Meteo_file_cell,METEO_FILE,'Delimiter',',')
            
            varargout.GENERATEMETEOFILEButton_3.Value = 0;
           varargout.newfilehasbeengeneratedLamp_3.Color = 'green';
            
        end

        if gen_qmelt_file_flag == 1
            
            varargout.newfilehasbeengeneratedLamp_4.Color = 'white';
            pause(0.1)
            
            %% QMELT
            if snowmelt_method == 0
                qmelt_estim_ts = max(TEMP_ts.data * T_index_coef,0);
                qmelt_estim_ts(isnan(qmelt_estim_ts)) = 0;
                qmelt_file = [time_pulse_sec,qmelt_estim_ts];
            else
               crhm_output = readtable(crhmoutput_dir);     
               timecrhm = str2double(crhm_output.time(2:end)) + 693960;
               snowmelt_data = str2double(crhm_output.snowmelt_int_1_); % m-> mm

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

            writecell(qmelt_file_cell,QMELT_FILE,'Delimiter',',')
            
           varargout.GENERATEQMELTFILEButton.Value = 0;
           varargout.newfilehasbeengeneratedLamp_4.Color = 'green';
            
        end

        %{
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
        %}

        try
            check = varargout.figure1.Position;
        catch
            clear all;
            winopen = 0;
        end

    end
end
