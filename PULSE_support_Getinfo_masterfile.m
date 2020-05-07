
% Retrieve data from masterfile
 function [comment,time_sim,H_LAY,L_LAY] = PULSE_support_Getinfo_masterfile(time_sim_elapsec,pulse_dir,masterfile)
    
    fid = fopen([pulse_dir,'/',masterfile]);
    
    while(~feof(fid))
        newline = fgetl(fid);
    
        if(contains(newline,'COMMNET'))
            comment = erase(newline,'COMMNET ');
        end
        if(contains(newline,'H_LAY'))
            H_LAY = str2double(erase(newline,'H_LAY '));
        end
         if(contains(newline,'L_LAY'))
            L_LAY = str2double(erase(newline,'L_LAY '));
        end
        if(contains(newline,'START_TIME'))
            start_time_text = erase(newline,'START_TIME ');
            start_time_num = datenum(start_time_text,'dd-mm-yyyy HH:MM:SS');
            time_sim = datenum(start_time_num + seconds(time_sim_elapsec));
        end  
    end
    
 end
 
