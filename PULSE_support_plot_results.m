
% Plot results
function PULSE_support_plot_results(pulse_dir,results_dir,chemical_species,...
                    col_li,masterfile,Obs_file)

    % load pulse results
     [time_sim_elapsec,h_layers_max,c_m,c_s,c_total,poros_m,poros_s,...
         v_liqwater,v_swe,v_air] = PULSE_support_load_pulse_results(results_dir,col_li);
                                        
    % get comment and timenum from masterfile
    [comment,time_sim,H_LAY,L_LAY] = PULSE_support_Getinfo_masterfile(time_sim_elapsec,pulse_dir,masterfile);
    %H_LAY = H_LAY/10; % mm to cm
    ih = 0:H_LAY:(h_layers_max-1)*H_LAY;
    
    % generate mesh grid for plot surf
    [Hmesh,Tmesh] = meshgrid(ih,time_sim);
    
    % Retrieve Obs data for c_total
    [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec,...
        cmax_ctotal] = PULSE_support_GetTrans_obs_data(Obs_file,c_total,chemical_species); 
   
    inanloc = isnan(c_total);
    
    Y_obs_mesh = Y_obs_mesh * 10; % cm to mm

    figure('name',comment)
    for i = 1:8
        
        if i==1; var_print = c_m; var_print_name = ['Concentration liquid phase (mg/l): ', chemical_species]; end
        if i==2; var_print = c_s; var_print_name = ['Concentration solid phase (mg/l): ', chemical_species]; end
        if i==3; var_print = c_total; var_print_name = ['Concentration snow (liquid + solid phases) (mg/l): ', chemical_species]; end
        
        if i==4; var_print = v_liqwater/(L_LAY); var_print_name = 'Volume liquid phase [mm/mm]'; end
        if i==5; var_print = v_swe/(L_LAY); var_print_name = 'Volume solid phase [mm/mm]'; end
        if i==6; var_print = v_air/(L_LAY); var_print_name = 'Volume air phase [mm/mm]'; end
        
        if i==7; var_print = poros_m; var_print_name = 'Volume fraction of liquid phase [-]'; end
        if i==8; var_print = poros_s; var_print_name = 'Volume fraction of solid phase [-]'; end
                
        var_print(inanloc) = NaN;
             
        subplot(3,3,i)
        surf(Tmesh,Hmesh,var_print)
        grid on
        if i == 3
            hold on
            h1 = scatter3(X_obs_mesh,...
                          Y_obs_mesh,...
                          Z_obs_mesh,...
                          Marsize_obs_mesh,...
                            colvec,...
                            'filled','markeredgecolor','b');
            cmax_i = cmax_ctotal;          
        %elseif i == 4 || i == 5
        %    cmax_i = 1;
        else
            cmax_i = max(max(var_print));
            cmin_i = min(min(var_print));
        end
        xlim([min(Tmesh(:,1)) max(Tmesh(:,1))])
        set(gca, 'XTick',linspace(min(Tmesh(:,1)),max(Tmesh(:,1)),6));
        set(gca, 'XTickLabel',linspace(min(Tmesh(:,1)),max(Tmesh(:,1)),6));
        ylim([0 max(Hmesh(1,:))])
        datetick('x','mmm-dd','keepticks','keeplimits')  
        try
            caxis([cmin_i cmax_i])
        catch
        end
        %colormap(othercolor('Blues9'))
        colormap(jet)
        colorbar
        view(0,90)
        %xlabel('Date')
        ylabel('Snow height [mm]')
        title(var_print_name,'Interpreter', 'none')
        shading interp
        set(gca, 'layer', 'top');
        alpha 0.7
        
    end
    
    Y_obs_mesh = Y_obs_mesh * 10; % cm -> mm
    
    %Model_data_interc = interp1(Data_time(1:end-1),Data(1:end-1),T_WQsort);
    Model_data_interc = interp2(Hmesh,Tmesh,c_total,Y_obs_mesh,X_obs_mesh);
    
    figure('name',[comment,' obsVSmodel'])
    scatter(Z_obs_mesh,Model_data_interc,'k')
    %limmin = min(min(Model_data_interc));
    limmax = max(max(Model_data_interc));
    hold on
    plot([0 limmax],[0 limmax],'k')
    xlim([0 limmax])
    ylim([0 limmax])
    grid on
    xlabel('Obs (mg/l)')
    ylabel('Model (mg/l)')
    
end

