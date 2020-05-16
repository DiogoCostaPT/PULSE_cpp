
% Retrieve data from masterfile
 function [comment,time_sim,H_LAY,L_LAY] = PULSE_support_Getinfo_masterfile(time_sim_elapsec,masterfile)
    
    fid = fopen(masterfile);
    
    while(~feof(fid))
        newline = fgetl(fid);
    
        if(contains(newline,'COMMNET'))
            comment = erase(newline,'COMMNET ');
        end
        if(contains(newline,'H_LAY_mm'))
            H_LAY = str2double(erase(newline,'H_LAY_mm '));
        end
         if(contains(newline,'L_LAY_mm'))
            L_LAY = str2double(erase(newline,'L_LAY_mm '));
        end
        if(contains(newline,'START_TIME'))
            start_time_text = erase(newline,'START_TIME ');
            start_time_num = datenum(start_time_text(2:end-1),'dd-mm-yyyy HH:MM:SS');
            time_sim = datenum(start_time_num + seconds(time_sim_elapsec));
        end  
    end
    
 end
 
