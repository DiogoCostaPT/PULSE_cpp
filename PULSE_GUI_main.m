

gui_varargout = PULSE_GUI();
winopen = 1;
pause(2);

while(winopen == 1)
    
    try
        GENALLFILES_FLAG = gui_varargout.GENERATEALLFILESButton.Value;
        genetate_masterfile = gui_varargout.GENERATEMASTERFILEButton_2.Value;
        gen_0txtfile_flag = gui_varargout.GENERATEICFILEButton_2.Value;   
        gen_meteo_file_flag = gui_varargout.GENERATEMETEOFILEButton_3.Value;
        gen_qmelt_file_flag = gui_varargout.GENERATEQMELTFILEButton.Value;
        RUN_single = gui_varargout.RUNButton.Value;
        RUN_sens = gui_varargout.RUNButton_2.Value;
    catch
        winopen = 0;
        continue;
    end
    
    execute_flag = GENALLFILES_FLAG...
                   + genetate_masterfile ...
                   + genetate_masterfile...
                   + gen_meteo_file_flag...
                   + gen_qmelt_file_flag...
                   + RUN_single...
                   + RUN_sens;
    
    if execute_flag == 0
        pause(2)
        continue; 
    end
    
    if GENALLFILES_FLAG == 1
        genetate_masterfile = 1;
        gen_0txtfile_flag = 1;
        gen_meteo_file_flag = 1;
        gen_qmelt_file_flag = 1;
    end
    
    % Master file
    masterfile_fullpath = gui_varargout.MASTERFILENAMEfullorrelativepathEditField.Value;
    COMMENT = gui_varargout.HYDRO_SOLVEREditField_2.Value;
    START_TIME = gui_varargout.START_TIMEEditField.Value;
    END_TIME = gui_varargout.END_TIMEEditField.Value;
    PRINT_STEP = gui_varargout.PRINT_STEPsecEditField.Value;
    L_LAY = gui_varargout.L_LAYmmEditField_2.Value; % mm
    H_LAY = gui_varargout.H_LAYmmEditField_2.Value; % 100 cm * 10 = 1000 mm
    VFRAC_AIR_FRESHSNOW = gui_varargout.VFRAC_AIR_FRESHSNOWEditField.Value;
    DENSITY_ICE = gui_varargout.DENSITY_ICEkgm3EditField.Value;
    DENSITY_WATER = gui_varargout.DENSITY_WATERkgm3EditField.Value;
    DENSITY_FRESHSNOW = gui_varargout.DENSITY_FRESHSNOWkgm3EditField.Value;
    A_D = gui_varargout.A_Dm2sEditField.Value;
    ALPHA_IE = gui_varargout.ALPHA_IEEditField.Value;
    COMPFACTOR = gui_varargout.COMPFACTOREditField.Value;
    HYDRO_SOLVER = gui_varargout.HYDRO_SOLVEREditField.Value;
    METEO_FILE = gui_varargout.METEO_FILEfullorrelativepathEditField.Value;
    QMELT_FILE = gui_varargout.QMELT_FILEfullorrelativepathEditField.Value;
    
    % IC file and metrofiles
    Results_folder_pulse = gui_varargout.SAVEFILEINFOLDERfullorrelativepathEditField.Value;
    T_field = gui_varargout.AirtemperatureEditField.Value;
    Rh_field = gui_varargout.RHEditField.Value;
    WS_field = gui_varargout.WindSpeedEditField.Value;
    QSI_field = gui_varargout.QsiEditField.Value;
    RAIN_field = gui_varargout.RainEditField.Value;
    SNOW_field = gui_varargout.SnowfallEditField.Value;
    crhmoutput_dir = gui_varargout.CRHMoutputfilefullorrelativepathEditField.Value;
    chemistry_file = gui_varargout.SnowChemistryfilefullorrelativepathEditField.Value;
    chemical_species = gui_varargout.ChemicalspeciesworksheetnameinthefileaboveEditField.Value; 
    H_SNOWPACK = gui_varargout.H_SNOWPACKmmEditField.Value; % 10 mm
    L_SNOWPACK = gui_varargout.L_SNOWPACKmmEditField.Value; % 10 mm
    v_frac_air_init = gui_varargout.VFRAC_AIR_FRESHSNOWEditField.Value; % volume percertage

    % qmelt and meteo files
    snowmelt_method = gui_varargout.snowmeltcalcmethod1Tindex2CRHMEditField_3.Value; % 0)T-index, 1) CRHM output
    T_index_coef = gui_varargout.T_index_coefonlyusedifsnowmeltcalcmethod1EditField.Value; % only used if snowmelt_method = 0
    %crhmoutput_dir = gui_varargout.CRHMoutputfilefullorrelativepathifsnowmeltcalcmethod2EditField.Value; % only used if snowmelt_method =  1;
    variableNameCRHM = gui_varargout.snowmeltevent_field.Value; % CRHM snowmelt variable name
    
    % reset uneditable boxes
    gui_varargout.SnowChemistryfilefullorrelativepathEditField_2.Value = gui_varargout.SnowChemistryfilefullorrelativepathEditField.Value;
    gui_varargout.ChemicalspeciesworksheetnameinthefileaboveEditField_2.Value = gui_varargout.ChemicalspeciesworksheetnameinthefileaboveEditField.Value; 
    gui_varargout.METEO_FILEfullorrelativepathEditField_2.Value = gui_varargout.METEO_FILEfullorrelativepathEditField.Value;
    gui_varargout.QMELT_FILEfullorrelativepathEditField_2.Value = gui_varargout.QMELT_FILEfullorrelativepathEditField.Value;
    gui_varargout.H_LAYmmEditField_2.Value = gui_varargout.H_LAYmmEditField.Value;
    gui_varargout.L_LAYmmEditField_2.Value = gui_varargout.L_LAYmmEditField.Value;
    gui_varargout.VFRAC_AIR_FRESHSNOWEditField_2.Value = gui_varargout.VFRAC_AIR_FRESHSNOWEditField.Value;
    
    % simuations
    pulse_dir = gui_varargout.Directoryofpulse_cppexecutableEditField.Value; % bin/
    col_li = round(L_SNOWPACK/L_LAY/2);

    Run_pulse_flag = 0;
    Sens_run_flag = 0;
    Clean_results_folder_except_IC_flag =0;
    Plot_results_flag = 0;
    Sens_analysis_flag = 0;
    Sens_plot_results_flag = 0;

    if RUN_single == 1
        Run_pulse_flag = gui_varargout.runsimulationCheckBox.Value;
        Clean_results_folder_except_IC_flag = gui_varargout.deleteResultsfolderexceptICfileCheckBox.Value;
        Plot_results_flag = gui_varargout.PlotresultssinglerunCheckBox.Value;
    elseif RUN_sens == 1
       Sens_run_flag = gui_varargout.RuntestsCheckBox.Value;
       Sens_analysis_flag = gui_varargout.AnalyzeresultsCheckBox.Value; %0;
       Sens_plot_results_flag = gui_varargout.PlotSensresultsCheckBox.Value;
    end

    %Run_pulse_flag = 1;
    IC_file = gui_varargout.edit3.Value; % 0.txt

    %Sens_run_flag = 0;
    num_samples = gui_varargout.scenariosEditField.Value; %500;
    A_D_max = gui_varargout.A_D_maxEditField.Value; %0.0003;           
    ALPHA_IE_max = gui_varargout.ALPHA_IE_maxEditField.Value; %0.000003;

    %% Run PULSE (once)
    if Run_pulse_flag
        RunningsimulationGauge.Value = 0;
        gui_varargout.Lamp_5.Color = [1 1 0];
        drawnow;
        
        PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,pulse_dir,...
        Results_folder_pulse,IC_file,masterfile_fullpath); 
        
        gui_varargout.RunningsimulationGauge.Value = 100;
        gui_varargout.Lamp_5.Color = 'green';
        drawnow;
    end

    %% Plot Results
    if Plot_results_flag
        gui_varargout.RunningsimulationGauge.Value = 0;
        gui_varargout.CalcsnowpackdepthGauge.Value = 0;
        gui_varargout.LoadingResultsGauge.Value = 0;
        gui_varargout.PlottingresultsGauge.Value = 0;
        gui_varargout.Lamp_6.Color = [1 1 0];
        gui_varargout.Lamp_7.Color = [1 1 0];
        gui_varargout.Lamp_8.Color = [1 1 0];
        drawnow;
        
        PULSE_support_plot_results(Results_folder_pulse,...
                chemical_species,col_li,masterfile_fullpath,chemistry_file); 
    
        gui_varargout.RunningsimulationGauge.Value = 100;
        gui_varargout.CalcsnowpackdepthGauge.Value = 100;
        gui_varargout.LoadingResultsGauge.Value = 100;
        gui_varargout.PlottingresultsGauge.Value = 100;
        gui_varargout.Lamp_6.Color = 'green';
        gui_varargout.Lamp_7.Color = 'green';
        gui_varargout.Lamp_8.Color = 'green';
        drawnow;
    end

    %% Sensitivity runs (multiple runs)
    if Sens_run_flag
        RunningscenariosGauge.Value = 0;
        gui_varargout.Lamp_2.Color = [1 1 0]; 
        drawnow;
        
        PULSE_support_Sens_analysys_run(masterfile_fullpath,pulse_dir,Results_folder_pulse,...
            num_samples,A_D_max,ALPHA_IE_max);
        
        gui_varargout.RunningscenariosGauge.Value = 100;
        gui_varargout.Lamp_2.Color = 'green';
        drawnow;
    end

    % Sensitivity analysis (processing the results in Sensitivity_analysis
    if Sens_analysis_flag
        gui_varargout.ProcessingresultsGauge.Value = 0;
        gui_varargout.Lamp_3.Color = [1 1 0];
        drawnow;
        
        PULSE_support_Sens_process(chemistry_file,chemical_species,col_li,...
            masterfile_fullpath);
        
        gui_varargout.ProcessingresultsGauge.Value = 100;
        gui_varargout.Lamp_3.Color = 'green';
        drawnow;
    end

    % Sensitivity results plotting
    if Sens_plot_results_flag
        gui_varargout.PlottingresultsGauge_2.Value = 0;
        gui_varargout.Lamp_4.Color = [1 1 0];
        drawnow;
        
        PULSE_support_Sens_plot(); 
        
        gui_varargout.PlottingresultsGauge_2.Value = 100;
        gui_varargout.Lamp_4.Color = 'green';
        drawnow;        
    end
    
    gui_varargout.RUNButton.Value = 0;
    gui_varargout.RUNButton_2.Value = 0;
    
    
    %% GENERATE MASTER FILE
    if genetate_masterfile == 1
        
        gui_varargout.newfilehasbeengeneratedLamp.Color = 'white';
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
        
        gui_varargout.GENERATEMASTERFILEButton_2.Value = 0;
        gui_varargout.newfilehasbeengeneratedLamp.Color = 'green';
        
    end
    
    %% Generate 0.txt file 
    if gen_0txtfile_flag == 1
        
        gui_varargout.newfilehasbeengeneratedLamp_2.Color = 'white';
        pause(0.1)

        dataraw_chem = xlsread(chemistry_file,chemical_species);

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
        
        folder0txt_path = [Results_folder_pulse,'/0.txt'];
        folder0txt_path = strrep(folder0txt_path,'//','/');

        writecell(file_0txt_cell,folder0txt_path,'Delimiter',',')
        
        gui_varargout.GENERATEICFILEButton_2.Value = 0;
        gui_varargout.newfilehasbeengeneratedLamp_2.Color = 'green';
        

    end
    
    
    %% Generate Qmelt data
    if gen_meteo_file_flag == 1 || gen_qmelt_file_flag == 1

        % chem data
        dataraw_chem = xlsread(chemistry_file,chemical_species);
        depths_obs = [10,20,30,40,50,60,70,80,90,100];
        time_obs_chem = dataraw_chem(:,1) + 693960;
        timestart_chem = time_obs_chem(1); % needs the start for initial conditions
        timeend_chem = datenum(END_TIME,'dd-mm-yyyy HH:MM:SS'); % end taken from simset
        
        % meteo data
        dataraw_meteo = readtable(crhmoutput_dir);
        dataraw_meteo = dataraw_meteo(2:end,:);
        time_meteo = str2double(dataraw_meteo.time) + 693960;

        % sub-set meteo data to chem data time period
        i_tstart = find(time_meteo==timestart_chem);
        i_tend = find(time_meteo==timeend_chem);
        dataraw_meteo_relsubset = dataraw_meteo(i_tstart:i_tend,:);

        % prepare data    
        TIME = str2double(dataraw_meteo_relsubset.time) + 693960;
        TIME_str = datestr(TIME);
        T_field = strrep(strrep(T_field,'(','_'),')','_');
        Rh_field = strrep(strrep(Rh_field,'(','_'),')','_');
        WS_field = strrep(strrep(WS_field,'(','_'),')','_');
        QSI_field = strrep(strrep(QSI_field,'(','_'),')','_');
        RAIN_field = strrep(strrep(RAIN_field,'(','_'),')','_');
        SNOW_field = strrep(strrep(SNOW_field,'(','_'),')','_');
        
        TEMP = str2double(dataraw_meteo_relsubset.(genvarname(T_field)));
        RH = str2double(dataraw_meteo_relsubset.(genvarname(Rh_field)));
        WS = str2double(dataraw_meteo_relsubset.(genvarname(WS_field)));
        RAD = str2double(dataraw_meteo_relsubset.(genvarname(QSI_field)));
        SNOWfall = str2double(dataraw_meteo_relsubset.(genvarname(SNOW_field)));
        RAIN = str2double(dataraw_meteo_relsubset.(genvarname(RAIN_field)));

        TEMP_ts = timeseries(TEMP,TIME);
        RH_ts = timeseries(RH,TIME);
        WS_ts = timeseries(WS,TIME);
        RAD_ts = timeseries(RAD,TIME);
        SNOWfall_ts = timeseries(SNOWfall,TIME);
        RAIN_ts = timeseries(RAIN,TIME);

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
        %RAIN_ts = PREC_ts;
        %RAIN_ts.Data(PREC_type_ts.Data~=2) = 0;

        %SNOWfall_ts = PREC_ts;
        %SNOWfall_ts.Data(PREC_type_ts.Data==2) = 0; % daily

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
            
            gui_varargout.newfilehasbeengeneratedLamp_3.Color = 'white';
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
            
            gui_varargout.GENERATEMETEOFILEButton_3.Value = 0;
           gui_varargout.newfilehasbeengeneratedLamp_3.Color = 'green';
            
        end

        if gen_qmelt_file_flag == 1
            
            gui_varargout.newfilehasbeengeneratedLamp_4.Color = 'white';
            pause(0.1)
            
            %% QMELT
            if snowmelt_method == 1
                qmelt_estim_ts = max(TEMP_ts.data * T_index_coef,0);
                qmelt_estim_ts(isnan(qmelt_estim_ts)) = 0;
                qmelt_file = [time_pulse_sec,qmelt_estim_ts];
            elseif snowmelt_method == 2
               crhm_output = readtable(crhmoutput_dir);     
               timecrhm = str2double(crhm_output.time(2:end)) + 693960;
               variableNameCRHM = strrep(variableNameCRHM,'(','_');
               variableNameCRHM = strrep(variableNameCRHM,')','_');
               snowmelt_data = str2double(crhm_output.(genvarname(variableNameCRHM))); % m-> mm

               date_start = datenum(START_TIME,'dd-mm-yyyy HH:MM:SS');
               date_end = timeend_chem;
               

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
            
           gui_varargout.GENERATEQMELTFILEButton.Value = 0;
           gui_varargout.newfilehasbeengeneratedLamp_4.Color = 'green';
            
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
        
       gui_varargout.GENERATEALLFILESButton.Value = 0;
        
        try
            check = gui_varargout.figure1.Position;
        catch
            clear all;
            winopen = 0;
        end

    end
end
