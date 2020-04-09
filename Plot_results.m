
% Plot Results

results_folder = '/media/dcosta/data/megasync/ec_main/models/pulse/code/code_pulse_cpp/Results/';

filenames_raw = dir(results_folder);

% get results
filenames = {filenames_raw.name};
filenames = filenames(contains(filenames,'.txt'));
filename_no = [];
for i=1:numel(filenames)
    file_i = filenames{i};
    time_i = str2double(file_i(1:end-4));
    filename_no = [filename_no,time_i];
end

filename_no_sort = sort(filename_no,'ascend');
filename_no_sort(isnan(filename_no_sort)) = [];

% decide the cell l to print
cell_li = 5;

% get the time steps
timesteps_num = numel(filename_no_sort);

% Get max height
hbar = parfor_progressbar(numel(filenames), 'Determine max height...');
h_layers = zeros(numel(filename_no_sort),1)*NaN;
parfor i = 1:numel(filename_no_sort)
    hbar.iterate(1)
    
    try
        dataraw = readtable([results_folder,num2str(filename_no_sort(i)),'.txt']);
        i_rel = find(dataraw.Var2==cell_li);
        data = dataraw(i_rel,:);
        h_layers(i) = numel(data(:,1));
    catch
        disp(['problem with results file ',num2str(filename_no_sort(i))])
    end
    
end
hbar.close();
h_layers_max = max(h_layers);

% create matrixes
time = ones(timesteps_num) * NaN;
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
        data_raw = readtable([results_folder,num2str(filename_no_sort(i)),'.txt']);
        
        time_i = str2double(file_i(1:end-4)); % in seconds
        
        i_rel = find(data_raw.Var2==cell_li);
        
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


figure
subplot(2,3,1)
surf(c_m)
ylim([0 timesteps_num])
xlim([0 h_layers_max])
colorbar
view(90,-90)
title('c_m')
shading interp
subplot(2,3,2)
surf(c_s)
colorbar
view(90,-90)
title('c_s')
shading interp
subplot(2,3,3)
surf(c_total)
colorbar
view(90,-90)
title('c_total')
shading interp
subplot(2,3,4)
surf(poros_m)
colorbar
view(90,-90)
title('poros_m')
shading interp
subplot(2,3,5)
surf(poros_s)
colorbar
view(90,-90)
title('poros_s')
shading interp