
% Run PULSE

function PULSE_support_run_pulse(Clean_results_folder_except_IC_flag,...
            pulse_dir,results_dir,IC_file,masterfile)

if Clean_results_folder_except_IC_flag
   deleted_results([results_dir,'/'],IC_file);
end

command = ['./',pulse_dir,'/pulse_cpp ',masterfile,' ',results_dir];
[status,cmdout] = system(command,'-echo');

end


% delete Results folder, except the 0.txt file
function deleted_results(results_dir,IC_file)
    f=dir(results_dir)
    f={f.name};
    n=find(strcmp(f,IC_file));
    f{n}=[];
    for k=1:numel(f)
      delete([results_dir,f{k}])
    end
end

