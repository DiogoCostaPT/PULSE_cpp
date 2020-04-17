
masterfile = 'simset.pulse';

Run_pulse_flag = 0;
Clean_results_folder_except_IC_flag = 0;
IC_file = '0.txt';

Plot_results = 1;
chemical_species = 'NO3';
col_li = 5; % vertical cell to print results
Obs_file = 'BRG_data.xlsx';

% directories
pulse_dir = 'bin';
results_dir = 'Results';


%% Run PULSE
if Run_pulse_flag
    
    if Clean_results_folder_except_IC_flag
       deleted_results('Results/',IC_file);
    end
    
    command = ['./',pulse_dir,'/pulse_cpp ',pulse_dir,'/',masterfile]
    [status,cmdout] = system(command,'-echo')
    
end

sec1 = datenum(now + seconds(1));
now - sec1

%% Plot Results
if Plot_results
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
    [time_sim_elapsec,c_m,c_s,c_total,poros_m,poros_s] = load_results(results_dir,filename_no_sort,...
                                            timesteps_num,h_layers_max,col_li);
                                        
    % get comment and timenum from masterfile
    [comment time_sim,H_LAY] = Getinfo_masterfile(time_sim_elapsec,pulse_dir,masterfile);
    H_LAY = H_LAY/10; % mm to cm
    ih = 0:H_LAY:(h_layers_max-1)*H_LAY;
    
    % generate mesh grid for plot surf
    [Hmesh,Tmesh] = meshgrid(ih,time_sim);
    
    % Retrieve Obs data for c_total
    [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec] = Get_obs_data(Obs_file,c_total,chemical_species); 
   

    figure('name',comment)
    for i = 1:5
        
        if i==1; var_print = c_m; var_print_name = 'c_m'; end
        if i==2; var_print = c_s; var_print_name = 'c_s'; end
        if i==3; var_print = c_total; var_print_name = 'c_total'; end
        if i==4; var_print = poros_m; var_print_name = 'poros_m'; end
        if i==5; var_print = poros_s; var_print_name = 'poros_s'; end
            
        subplot(2,3,i)
        surf(Tmesh,Hmesh,var_print)
        if i == 3
            hold on
            h1 = scatter3(X_obs_mesh,...
                          Y_obs_mesh,...
                          Z_obs_mesh,...
                          Marsize_obs_mesh,...
                            colvec,...
                            'filled','markeredgecolor','w');
            cmax_i = cmax_ctotal;          
        else
            cmax_i = max(max(var_print));
        end
        xlim([min(Tmesh(:,1)) max(Tmesh(:,1))])
        datetick('x','mm-dd','keepticks','keeplimits')
        ylim([0 max(Hmesh(1,:))])
        caxis([0 cmax_i])
        colormap(jet)
        colorbar
        view(0,90)
        title(var_print_name)
        shading interp
        
    end
    
end

%% AUXILIARY FUNCTIONS

% get num value
function num_i = getnums(str_i)

    num_i = [];
    for l = 1:numel(str_i)
        num_i = [num_i, str2double(regexp(str_i{l},'\d*','match','once'))];
    end
    
end

% delete Results folder, except the 0.txt file
function deleted_results(results_dir,IC_file)
    f=dir(results_dir)
    f={f.name}
    n=find(strcmp(f,IC_file));
    f{n}=[]
    for k=1:numel(f)
      delete([results_dir,'/',f{k}])
    end
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
            i_rel = dataraw.Var2==col_li;
            data = dataraw(i_rel,:);
            h_layers(i) = numel(data(:,1));
        catch
            disp(['Error: problem reading file ',file_2_read])
        end

    end
    h_layers_max = max(h_layers);
    hbar.close();
 end
 
 % Load results
 function [time,c_m,c_s,c_total,poros_m,poros_s] = load_results(results_dir,filename_no_sort,...
                                                    timesteps_num,h_layers_max,col_li)
    
    % create matrixes
    time = ones(timesteps_num,1) * NaN;
    ih = ones(timesteps_num,1) * NaN;
    c_m = ones(timesteps_num,h_layers_max) * NaN;
    c_s = ones(timesteps_num,h_layers_max) * NaN;
    poros_m = ones(timesteps_num,h_layers_max) * NaN;
    poros_s = ones(timesteps_num,h_layers_max) * NaN;

    % load data and process it
    hbar = parfor_progressbar(numel(filename_no_sort), 'Loading results...');
    for i = 1:numel(filename_no_sort)
        hbar.iterate(1)
        file_i = filename_no_sort(i);
        try
            data_raw = readtable([results_dir,'/',num2str(filename_no_sort(i)),'.txt']);

            time_i = file_i; % in seconds

            i_rel = find(data_raw.Var2==col_li);

            data = data_raw(i_rel,:);
            h_layers_i = numel(data(:,1));

            time(i) = time_i;
            
            c_m(i,1:h_layers_i) = flipud(data.Var3)';
            c_s(i,1:h_layers_i) = flipud(data.Var5)';
            poros_m(i,1:h_layers_i) = flipud(data.Var6)';
            poros_s(i,1:h_layers_i) = flipud(data.Var7)';

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
 
 % Retrieve data from masterfile
 function [comment,time_sim,H_LAY] = Getinfo_masterfile(time_sim_elapsec,pulse_dir,masterfile)
    
    fid = fopen([pulse_dir,'/',masterfile]);
    
    while(~feof(fid))
        newline = fgetl(fid);
    
        if(contains(newline,'COMMNET'))
            comment = erase(newline,'COMMNET ');
        end
        if(contains(newline,'H_LAY'))
            H_LAY = str2double(erase(newline,'H_LAY '));
        end
        if(contains(newline,'START_TIME'))
            start_time_text = erase(newline,'START_TIME ');
            start_time_num = datenum(start_time_text,'dd-mm-yyyy HH:MM:SS');
            time_sim = datenum(start_time_num + seconds(time_sim_elapsec));
        end  
    end
    
 end
 
 % Get obs data
 function [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec] = Get_obs_data(Obs_file,c_total,chemical_species)
   
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
    [X,Y] = meshgrid(depth_fixed_int,time);
    Z = data;
    C = data;

    depth_corr_onecolmn = reshape(depth_corr,[],1);
    Y_onecolmn = reshape(Y,[],1);
    Z_onecolmn = reshape(Z,[],1);
    C_onecolmn = reshape(C,[],1);
    
    X_obs_mesh = repmat(time,numel(elev_meas),1);
    Y_obs_mesh = depth_corr_onecolmn;
    Z_obs_mesh = Z_onecolmn + 20;
    Marsize_obs_mesh = 30*ones(numel(Z_onecolmn),1);
    
    cmax_obs = max(C_onecolmn);
    cmax_mod = max(max(c_total));
    cmax_ctotal = max(cmax_obs,cmax_mod);
    hotbar = jet;
    %hotbar = [hotbar];
    C_onecolmn(isnan(C_onecolmn)) = 0;
    colorc = round(C_onecolmn/cmax_ctotal * (numel(hotbar(:,1))-1))+1;
    colvec = hotbar(colorc,:);
    
 end