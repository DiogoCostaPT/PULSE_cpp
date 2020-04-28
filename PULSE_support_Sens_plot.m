
function PULSE_support_Sens_plot()

sens_folder = [pwd,'/Sensitivity_analysis/'];
windowtitle = 'Select the folder with the Sensitivity batch runs of interest'; 
sensbatch_dir = uigetdir(sens_folder,windowtitle);


sens_res_fullpath = [sensbatch_dir,'/Sens_results.mat'];
load(sens_res_fullpath)


figure
subplot(1,3,1)
ilox_best = find(Nash == max(Nash));
scatter3(A_D_all,ALPHA_IE_all,Nash,'k','marker','.')
hold on
scatter3(A_D_all(ilox_best),ALPHA_IE_all(ilox_best),Nash(ilox_best),'ro','filled')
xlabel('A_D (x)')
ylabel('ALPHA_IE (y)')
zlabel('Nash (z)')
%view(0,90)
subplot(1,3,2)
ilox_best = find(RMSE == min(RMSE));
scatter3(A_D_all,ALPHA_IE_all,RMSE,'k','marker','.')
hold on
scatter3(A_D_all(ilox_best),ALPHA_IE_all(ilox_best),RMSE(ilox_best),'ro','filled')
xlabel('A_D (x)')
ylabel('ALPHA_IE (y)')
zlabel('RMSE (z)')
%view(0,90)
subplot(1,3,3)
ilox_best = find(Bias == min(Bias));
scatter3(A_D_all,ALPHA_IE_all,Bias,'k','marker','.')
hold on
scatter3(A_D_all(ilox_best),ALPHA_IE_all(ilox_best),Bias(ilox_best),'ro','filled')
xlabel('A_D (x)')
ylabel('ALPHA_IE (y)')
zlabel('Bias (z)')
%view(0,90)

end