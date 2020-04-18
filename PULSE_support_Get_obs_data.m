% get obs filedata
function [depth_fixed_int,depth_corr,time,data,elev_meas] = PULSE_support_Get_obs_data(Obs_file,chemical_species)

    BRG_data = importdata(Obs_file);
    datafields = fieldnames(BRG_data.data);
    depth = BRG_data.data.(genvarname(datafields{1}))(:,2:end-1); % last column is the accumulation average

    depth_temporar_increment = 20; % needed to make sure that the interpolation method will work
    depth_corr = -depth + max(100);

    depth_corr_max = max(max(depth_corr));
    depth_corr_min = min(min(depth_corr));

    depth_fixed_int = [depth_corr_min:0.1:depth_corr_max];

    iloc = find(strcmp(datafields,chemical_species));

    depths_str_baseline = BRG_data.textdata.(genvarname(datafields{iloc}));
    depths_raw = getnums(depths_str_baseline);   
    depths_cols_noNaN = find(~isnan(depths_raw)==1);


    % Get matrix data    
    data_raw = BRG_data.data.(genvarname(datafields{iloc}));
    time = data_raw(:,1) + 693960;
    data = data_raw(:,depths_cols_noNaN+1)/1000; %ppb to mg/l
    depths = depths_raw(depths_cols_noNaN);
    elev_meas = max(depths) - depths;
    
end

 
 % get num value
function num_i = getnums(str_i)

    num_i = [];
    for l = 1:numel(str_i)
        num_i = [num_i, str2double(regexp(str_i{l},'\d*','match','once'))];
    end
    
end
