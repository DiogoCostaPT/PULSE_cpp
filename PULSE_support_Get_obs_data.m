% get obs filedata
function [depth_fixed_int,depth_corr,time,data,elev_meas] = PULSE_support_Get_obs_data(Obs_file,chemical_species)

    %BRG_data = readtable(Obs_file,'Sheet','depth_correction');
    %datafields = fieldnames(BRG_data.data);
    depth = readtable(Obs_file,'Sheet','depth_correction'); % last column is the accumulation average
    %depth = table2array(depth_timetable);
    depth = depth(:,1:end-1);
    
    depth_temporar_increment = 20; % needed to make sure that the interpolation method will work
    depth_corr = - table2array(depth(:,2:end)) + 100;

    depth_corr_max = max(max(depth_corr));
    depth_corr_min = min(min(depth_corr));

    depth_fixed_int = [depth_corr_min:0.1:depth_corr_max];

    %iloc = find(strcmp(datafields,chemical_species));

    %depths_str_baseline = BRG_data.textdata.(genvarname(datafields{iloc}));
    depths_raw = [10,20,30,40,50,60,70,80,90,100];   
    depths_cols_noNaN = find(~isnan(depths_raw)==1);


    % Get matrix data    
    data_raw = readtable(Obs_file,'Sheet',chemical_species);
    time = datenum(data_raw.Var1);
    data = table2array(data_raw(:,2:end))/1000; %ppb to mg/l
    depths = depths_raw;
    elev_meas = max(depths) - depths;
    
end

 
 % get num value
function num_i = getnums(str_i)

    num_i = [];
    for l = 1:numel(str_i)
        num_i = [num_i, str2double(regexp(str_i{l},'\d*','match','once'))];
    end
    
end
