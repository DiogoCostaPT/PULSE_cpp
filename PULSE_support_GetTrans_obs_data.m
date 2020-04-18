
 % Get obs data
 function [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec,...
     cmax_ctotal] = PULSE_support_GetTrans_obs_data(Obs_file,c_total,chemical_species)
   
    [depth_fixed_int,depth_corr,time,data,elev_meas] = PULSE_support_Get_obs_data(Obs_file,chemical_species);    
    
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

