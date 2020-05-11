

function PULSE_support_Sens_analysys_run(masterfile,pulse_dir,num_samples,...
        A_D_max,ALPHA_IE_max)
    
sensfoldername = 'Sensitivity_analysis';
sens_dir = [pwd,'/',sensfoldername];

% create "Sensitivity_analysis" dir if inexistent
folders_raw = dir(pwd);
folders = {folders_raw.name};
isSensdir = find(strcmp(folders,sensfoldername));

if isempty(isSensdir)
    mkdir(sens_dir);
    test_num_max = 0;
else
    cd(sens_dir)
    sensdir_foldnames_raw = dir(sens_dir);
    cd ..
    sensdir_foldnames = {sensdir_foldnames_raw.name};
    test_num_all = [];
    for i=1:numel(sensdir_foldnames)
        try
            test_foldname = sensdir_foldnames{i};
            test_num = str2double(test_foldname(6:end));
            test_num_all = [test_num_all,test_num];
        catch
        end
    end
    test_num_all(isnan(test_num_all)) = [];
    test_num_max = max(test_num_all);
    if isempty(test_num_max); test_num_max = 0; end
end

new_test_dir = [sens_dir,'/Sens_',num2str(test_num_max+1)];
mkdir(new_test_dir);

A_D_all = rand(num_samples,1) * A_D_max;           
ALPHA_IE_all = rand(num_samples,1) * ALPHA_IE_max;

hbar = parfor_progressbar(num_samples, 'Sensitivity_test running...');
    
parfor i=1:num_samples
    
    A_D_all_i = A_D_all(i);
    ALPHA_IE_all_i = ALPHA_IE_all(i);
    
    % create run folder
    sim_subfold = [new_test_dir,'/Run_',num2str(i)];
    mkdir(sim_subfold);
        
    % copy pulse_cpp and input files
    mastfile_fullpath = [pulse_dir,'/',masterfile];
    copyfile(mastfile_fullpath,sim_subfold)
    pulse_fullpath = [pulse_dir,'/pulse_cpp'];
    copyfile(pulse_fullpath,sim_subfold)
    
    % extract input files from masterfile,copy to test folder, update
    % master file with new parameter and removing subfolders
    [qmelt_fullpath,meteo_fullpath] = get_inputfiles(mastfile_fullpath);
    copyfile(qmelt_fullpath,sim_subfold)
    copyfile(meteo_fullpath,sim_subfold)
    update_masterfile(sim_subfold,masterfile,A_D_all_i,ALPHA_IE_all_i);
    
    % create Results folder
    res_subfold = [sim_subfold,'/Results'];
    mkdir(res_subfold);

    % copy 0.txt file (IC conditions)
    file_0txt_fullpath = [pulse_dir,'/Results/0.txt'];
    copyfile(file_0txt_fullpath,res_subfold)
    
    % save sensitivity analysis info
    sens_info = {['A_D ', num2str(A_D_all_i)];
                ['ALPHA_IE ', num2str(ALPHA_IE_all_i)]};
    writecell(sens_info,[sim_subfold,'/Sens_info']);
    
    % run scenario
    PULSE_support_run_pulse(0,sim_subfold(numel(pwd)+2:end),...
        res_subfold,[],masterfile)
    
    hbar.iterate(1);   % update progress by one iteration
end
close(hbar);   %close progress bar


end

% get input files from master file
function [QMELT_FILE,METEO_FILE] = get_inputfiles(mastfile_fullpath)

fid = fopen(mastfile_fullpath);

 while(~feof(fid))
        newline = fgetl(fid);
            
        if(contains(newline,'QMELT_FILE'))
            QMELT_FILE = erase(newline,'QMELT_FILE ');
        end
        if(contains(newline,'METEO_FILE'))
            METEO_FILE = erase(newline,'METEO_FILE ');
        end
 end
fclose(fid)

end

% update master file (remove any subfolder in the path)
function update_masterfile(sim_subfold,masterfile,A_D_all_i,ALPHA_IE_all_i)

new_masterfile_dir = [sim_subfold,'/',masterfile];
fid = fopen(new_masterfile_dir);
new_masterfile = {};

 while(~feof(fid))
        newline = fgetl(fid);
            
        if(contains(newline,'QMELT_FILE'))
            newline = erase(newline,'QMELT_FILE ');
            iloc = strfind(newline,'/')+1;
            newline = newline(iloc:end);
            newline = ['QMELT_FILE ',sim_subfold,'/',newline];
        end
        if(contains(newline,'METEO_FILE'))
            newline = erase(newline,'METEO_FILE ');
            iloc = strfind(newline,'/')+1;
            newline = newline(iloc:end);
            newline = ['METEO_FILE ',sim_subfold,'/',newline];
        end
        if(contains(newline,'A_D'))
            newline = ['A_D ',num2str(A_D_all_i)];
        end
        if(contains(newline,'ALPHA_IE'))
            newline = ['ALPHA_IE ',num2str(ALPHA_IE_all_i)];
        end
        
        new_masterfile=[new_masterfile;newline];
        
end
fclose(fid);


fid = fopen(new_masterfile_dir,'W');
for i=1:numel(new_masterfile)
    new_masterfile_line = new_masterfile{i};
    formatSpec = ['%',num2str(numel(new_masterfile_line)),'s\n'];
    fprintf(fid,formatSpec,new_masterfile_line);
end
fclose(fid);


end


