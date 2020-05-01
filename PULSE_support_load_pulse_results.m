
function [time_sim_elapsec,h_layers_max,c_m,c_s,c_total,poros_m,poros_s,...
    v_liqwater,v_swe,v_air] = PULSE_support_load_pulse_results(results_dir,col_li)

filenames_raw = dir(results_dir);

    % get results' timesteps
    filenames = {filenames_raw.name};
    filenames = filenames(contains(filenames,'.txt'));
    filename_no = [];
    for i=1:numel(filenames)
        file_i = filenames{i};
        time_i = str2double(file_i(1:end-4));
        filename_no = [filename_no,time_i];
    end
    
    % sort results
    filename_no_sort = sort(filename_no,'ascend');
    filename_no_sort(isnan(filename_no_sort)) = [];

    % get the num of time steps
    timesteps_num = numel(filename_no_sort);

    % Get max snow depth
    h_layers_max = Get_max_snow_depth(results_dir,filename_no_sort,col_li);
      
    % load results
    [time_sim_elapsec,c_m,c_s,c_total,poros_m,poros_s,...
        v_liqwater,v_swe,v_air] = read_results(results_dir,filename_no_sort,...
                                            timesteps_num,h_layers_max,col_li);
                                        
end

% Get max snow depth
 function h_layers_max = Get_max_snow_depth(results_dir,filename_no_sort,col_li)
    hbar = parfor_progressbar(numel(filename_no_sort), 'Determine max snowpack depth...');
    h_layers = zeros(numel(filename_no_sort),1)*NaN;
    parfor i = 1:numel(filename_no_sort)
        hbar.iterate(1)

        try
            file_2_read = [results_dir,'/',num2str(filename_no_sort(i)),'.txt'];
            dataraw = readtable(file_2_read);
            if ~isempty(dataraw)
                i_rel = dataraw.Var2==col_li;
                data = dataraw(i_rel,:);
                h_layers(i) = numel(data(:,1));
            else
                h_layers(i) = 0;
            end
        catch
            disp(['Error: problem reading file ',file_2_read])
        end

    end
    h_layers_max = max(h_layers);
    hbar.close();
 end


 % Load results
 function [time,c_m,c_s,c_total,poros_m,poros_s,v_liqwater,v_swe,v_air] = read_results(results_dir,filename_no_sort,...
                                                    timesteps_num,h_layers_max,col_li)
    
    % create matrixes
    time = ones(timesteps_num,1) * NaN;
    %ih = ones(timesteps_num,1) * NaN;
    c_m = ones(timesteps_num,h_layers_max) * NaN;
    c_s = ones(timesteps_num,h_layers_max) * NaN;
    poros_m = ones(timesteps_num,h_layers_max) * NaN;
    poros_s = ones(timesteps_num,h_layers_max) * NaN;
    v_liqwater = ones(timesteps_num,h_layers_max) * NaN;
    v_swe = ones(timesteps_num,h_layers_max) * NaN;
    v_air = ones(timesteps_num,h_layers_max) * NaN;
    
    nh_l = numel(c_m(1,:));

    % load data and process it
    hbar = parfor_progressbar(numel(filename_no_sort), 'Loading results...');
    for i = 1:numel(filename_no_sort)
        hbar.iterate(1)
        file_i = filename_no_sort(i);
        try
            data_raw = readtable([results_dir,'/',num2str(filename_no_sort(i)),'.txt']);

            time_i = file_i; % in seconds
            
            if ~isempty(data_raw)
                i_rel = find(data_raw.Var2==col_li);

                data = data_raw(i_rel,:);
                h_layers_i = numel(data(:,1));

                time(i) = time_i;

                cm_i = flipud(data.Var3)';
                nh_i = numel(cm_i);
                extra_h = nh_l-nh_i;
                c_m(i,:) = [cm_i,zeros(1,extra_h)];

                cs_i = flipud(data.Var5)'; 
                c_s(i,:) =  [cs_i,zeros(1,extra_h)];

                poros_m_i = flipud(data.Var6)'; 
                poros_m(i,:) =  [poros_m_i,zeros(1,extra_h)];

                poros_s_i = flipud(data.Var7)'; 
                poros_s(i,:) =  [poros_s_i,zeros(1,extra_h)];
                
                v_liqwater_i = flipud(data.Var8)'; 
                v_liqwater(i,:) = [v_liqwater_i,zeros(1,extra_h)];
                
                v_swe_i = flipud(data.Var9)'; 
                v_swe(i,:) = [v_swe_i,zeros(1,extra_h)];
                
                v_air_i = flipud(data.Var10)'; 
                v_air(i,:) = [v_air_i,zeros(1,extra_h)];
                
            else
                 c_m(i,:) = NaN;
                 c_s(i,:) = NaN;
                 poros_m(i,:) = NaN;
                 poros_s(i,:) = NaN;
                 v_liqwater(i,:) = NaN;
                 v_swe(i,:) = NaN;
                 v_air(i,:) = NaN;
                 
            end

        catch
            disp(['problem with results file ',num2str(file_i)])
        end


    end
    close(hbar)

    c_total = c_s .* poros_s + c_m .* poros_m;

    c_m(c_m==0) = NaN;
    c_s(c_s==0) = NaN;
    c_total(c_total==0) = NaN;
    
    
 end