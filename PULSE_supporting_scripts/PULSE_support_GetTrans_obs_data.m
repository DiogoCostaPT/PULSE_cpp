% Andrea data
 function [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh] = PULSE_support_GetTrans_obs_data(Obs_file,chemical_species)
   
    [depth_corr,time,data,elev_meas] = PULSE_support_Get_obs_data(Obs_file,chemical_species);    
    
    %[X,Y] = meshgrid(depth_fixed_int,time);
    Z = data;
    C = data;

    depth_corr_onecolmn = reshape(depth_corr,[],1);
    %Y_onecolmn = reshape(Y,[],1);
    Z_onecolmn = reshape(Z,[],1);
    C_onecolmn = reshape(C,[],1);
    
    X_obs_mesh = repmat(time,numel(elev_meas),1);
    Y_obs_mesh = depth_corr_onecolmn;
    Z_obs_mesh = Z_onecolmn;
    Marsize_obs_mesh = 30*ones(numel(Z_onecolmn),1);
    
    
 end