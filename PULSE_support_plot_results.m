
% Plot results
function PULSE_support_plot_results(results_dir,chemical_species,...
                    col_li,masterfile,Obs_file)

    Obs_WQ_file_smooth = '/media/dcosta/data/megasync/ec_main/models/pulse/code/code_snowpack/snowpack_pulse/tools/BRG_data_smoothed.mat';
    smooth_mat_data_flag = 1;            
                
    obs_conc_force_max_caxis = {'true',1, 6, 1, 15};
    chem_i = 4;
    num_bins = 200;
    
    
    % load pulse results
     [time_sim_elapsec,h_layers_max,c_m,c_s,c_total,poros_m,poros_s,...
         v_liq,v_swe,v_air] = PULSE_support_load_pulse_results(results_dir,col_li);
                                        
    % get comment and timenum from masterfile
    [comment,time_sim,H_LAY,L_LAY] = PULSE_support_Getinfo_masterfile(time_sim_elapsec,masterfile);
    %H_LAY = H_LAY/10; % mm to cm
    ih = 0:H_LAY:(h_layers_max-1)*H_LAY;
    
    % generate mesh grid for plot surf
    [Hmesh,Tmesh] = meshgrid(ih,time_sim);
    
    % Retrieve Obs data for c_total
    Obs_exist_flag = false;
    
    %{
    try
        [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh,colvec,...
            cmax_ctotal] = PULSE_support_GetTrans_obs_data(Obs_file,c_total,chemical_species); 
        %Y_obs_mesh = Y_obs_mesh * 10; % cm to mm
        Obs_exist_flag = true;
    catch
    end
    %}
    
    if ~smooth_mat_data_flag
        [X_obs_mesh,Y_obs_mesh,Z_obs_mesh,Marsize_obs_mesh] = PULSE_support_GetTrans_obs_data(Obs_WQ_file,chemical_species);
                
     else % using smoothed data
        load(Obs_WQ_file_smooth)
        X_obs_mesh = obs_time_extract.(genvarname(chemical_species));
        Y_obs_mesh = obs_depth_extract.(genvarname(chemical_species));
        Z_obs_mesh = obs_data_extract.(genvarname(chemical_species));
        Marsize_obs_mesh = 30 * ones(size(obs_data_extract));
     end

    Z_obs_mesh(isnan(Z_obs_mesh)) = 0;
    Z_obs_mesh( Z_obs_mesh<0 ) = 0;

    if ~strcmp(obs_conc_force_max_caxis{1},'true')
        cmax_obs = max(max(Z_obs_mesh));
        cmax_mod = max(max(Z));
        cmax_ctotal = max(cmax_obs,cmax_mod);
    else
        cmax_ctotal = obs_conc_force_max_caxis{chem_i+1};
    end
    
    colorcode = 'YlOrRd'; % 'Blues' - liquid content, 'RdBu' - temperature
    colorm = brewermap(num_bins,colorcode);
    
    colorc = round(Z_obs_mesh/cmax_ctotal * (numel(colorm(:,1))-1))+1;
    for e = 1:numel(colorc)
        colorc(e) = min(num_bins, colorc(e));
    end
    colvec = colorm(colorc,:);

    
   
    inanloc = isnan(c_total);
    
    Hmesh = Hmesh/10; % mm -> cm

    figure('name',comment)
    for i = 1:8
        
        if i==1; var_print = c_m; var_print_name = ['Concentration liquid phase (mg/L): ', chemical_species]; end
        if i==2; var_print = c_s; var_print_name = ['Concentration solid phase (mg/L): ', chemical_species]; end
        if i==3; var_print = c_total; var_print_name = ['Concentration snow (liquid + solid phases) (mg/L): ', chemical_species]; end
        
        if i==4; var_print = v_liq; var_print_name = 'Volume liquid phase [mm/mm/m]'; end
        if i==5; var_print = v_swe; var_print_name = 'Volume solid phase [mm/mm/m]'; end
        if i==6; var_print = v_air; var_print_name = 'Volume air phase [mm/mm/m]'; end
        
        if i==7; var_print = poros_m; var_print_name = 'Volume fraction of liquid phase [-]'; end
        if i==8; var_print = poros_s; var_print_name = 'Volume fraction of solid phase [-]'; end
                
        var_print(inanloc) = NaN;
        
        var_print(var_print<0) = 0;
             
        subplot(3,3,i)
        
        surf(Hmesh,Tmesh,var_print)
        grid on
        if i == 3
            hold on
            %if Obs_exist_flag
            Y_obs_mesh = 111 - Y_obs_mesh;
                h1 = scatter3(Y_obs_mesh,...
                              X_obs_mesh,...
                              -1*ones(size(Y_obs_mesh)),...
                              Marsize_obs_mesh,...
                                colvec,...
                                'filled','markeredgecolor','k');
            %end
            cmax_i = cmax_ctotal;          
        %elseif i == 4 || i == 5
        %    cmax_i = 1;
        else
            cmax_i = max(max(var_print));
            cmin_i = min(min(var_print));
        end
        ylim([min(Tmesh(:,1)) max(Tmesh(:,1))])
        set(gca, 'YTick',linspace(min(Tmesh(:,1)),max(Tmesh(:,1)),6));
        set(gca, 'YTickLabel',linspace(min(Tmesh(:,1)),max(Tmesh(:,1)),6));
        xlim([0 max(Hmesh(1,:))])
        datetick('y','mmm-dd','keepticks','keeplimits')  
        try
            caxis([cmin_i cmax_i])
        catch
        end
        %colormap(othercolor('Blues9'))
        colorcode = 'YlOrRd'; % 'Blues' - liquid content, 'RdBu' - temperature
        colorm = brewermap(num_bins,colorcode);
        %colorm(1,1) = 0.8;
        %colorm(1,2) = 0.8;
        %colorm(1,3) = 0.8;
        colormap(colorm)
        hcb = colorbar;
        
        view(90,-90)
        %xlabel('Date')
        xlabel('Snowdepth [cm]')
        title(var_print_name,'Interpreter', 'none')
        shading interp
        set(gca, 'layer', 'top');
        alpha 0.7
        
        
    end
    
    Y_obs_mesh = Y_obs_mesh; % cm -> mm
    
    %Model_data_interc = interp1(Data_time(1:end-1),Data(1:end-1),T_WQsort);
    Model_data_interc = interp2(Hmesh,Tmesh,c_total,Y_obs_mesh,X_obs_mesh);
    
    %[r2 rmse] = rsquare(Z_obs_mesh,Model_data_interc);
    
     % remove NaNs
     Model_data_interc_use = Model_data_interc;
     WQ_use = Z_obs_mesh;

     
    nanloc = find(isnan(Model_data_interc_use));
    Model_data_interc_use(nanloc) = [];
    WQ_use(nanloc) = [];
    
    nanloc = find(isnan(WQ_use));
    Model_data_interc_use(nanloc) = [];
    WQ_use(nanloc) = [];
    
    %figure
    %scatter(Z_obs_mesh,Model_data_interc)
    %limmin = min(min(Model_data_interc));
    %limmax = max(max(Model_data_interc));
    %hold on
    %plot([limmin limmin],[limmax limmax],'k')
    %xlim([limmin limmax])
    %grid on
    %xlabel('Obs (mg/l)')
    %ylabel('Model (mg/l)')
    
    
    % Nash
    numerator=(WQ_use-Model_data_interc_use).^2;
    denominator=(WQ_use-mean(WQ_use)).^2;
    Nash =1-(sum(numerator)/sum(denominator));
    
    % RMSE
    Sumcal = (Model_data_interc_use-WQ_use).^2;
    numerator = sum(Sumcal);
    n=numel(WQ_use);
    RMSE=(numerator/n)^(1/2);
    
    % BIAS
    numerator = sum(WQ_use);
    denominator = sum(Model_data_interc_use);
    Bias = numerator/denominator-1;
    
    
    figure('name',[comment,' obsVSmodel'])
    scatter(Z_obs_mesh,Model_data_interc,'k')
    %limmin = min(min(Model_data_interc));
    limmax = max(max(Model_data_interc));
    hold on
    plot([0 limmax],[0 limmax],'k')
    %xlim([0 limmax])
    %ylim([0 limmax])
    axis tight
    grid on
    xlabel('Obs (mg/l)')
    ylabel('Model (mg/l)')
    title(['RMSE = ', num2str(RMSE),'; Nash = ',num2str(Nash),'; Bias = ',num2str(Bias)])
    
    % Export
    export_ts = timeseries(c_m(:,2),datestr(time_sim));
    figure
    plot(export_ts,'r-','linewidth',1)
    grid on
    title({'Meltwater discharge concentration: ', chemical_species})
    ylabel('mg/L')
    
    
end

function [r2 rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).
if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end
% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end
% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);
if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end
rmse = sqrt(mean((y(:) - f(:)).^2));
end