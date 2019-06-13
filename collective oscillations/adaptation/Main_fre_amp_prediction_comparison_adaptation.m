running_fresh_computation=1;
if running_fresh_computation
    
    disp("running_fresh_computation!!!")
    plot_on=1;
    
    
% which_model=1; % nonlinear linear signal + linear (negative-feedback) adaptation
 which_model=2; % nonlinear negative-feedback adaptation + linear signal
% which_model=3; % nonlinear feedforward adaptation + linear signal


%% parameters
tau_s=1;
tau_a=1;
tau_y=1;
gamma_a=1;
epsilon=0.1;
alpha1=0.5;
alpha2=0.5;

para_array=[tau_s,tau_a,tau_y,gamma_a,epsilon,alpha1,alpha2];

%% oscillation amplitude and frequency from directly simulating the collective dynamics
N_min=floor(2/(alpha1*alpha2));
N_array=N_min:1:(N_min+5); % range of the cell density used in the paper, should be integer
[N_bar,os_fre_h,A_array_simu,B_array_simu]=computing_all_oscillation_amp_fre_adaptation(N_array,para_array,which_model,plot_on);


%% predicted oscillation amplitude and frequency,  
Pert_array=0:0.1:3; %regime of amplitude perturbation
[predicted_N,predicted_A,predicted_B,predicted_w]=predicting_amp_fre_from_analytical_response_spectrum_adaptation(which_model,para_array,Pert_array,plot_on);


%% save_all_results

% tt=today;
% aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/amp_fre_predicted_comparison_adaptation";
save(output_name);

else
    load ./data/amp_fre_predicted_comparison_adaptation

end



figureParameter
f1=plot(N_bar,A_array_simu,'*r',predicted_N,predicted_A,'-b');
xlim([0 5]);
%ylim([-1 1]);
a1=xlabel('$\bar{N}$');
a2=ylabel('$A$');
h1=legend("Simulation","Theory")
set(h1,'location','southwest')
fig_name='./figure/N_A.eps';
figurePostTreat


figureParameter
f1=plot(N_bar,B_array_simu,'*r',predicted_N,predicted_B,'-b');
xlim([0 5]);
%ylim([0 20]);
a1=xlabel('$\bar{N}$');
a2=ylabel('$B$');
h1=legend("Simulation","Theory")
set(h1,'location','southwest')
fig_name='./figure/N_B.eps';
figurePostTreat

% figureParameter
% f1=plot(N_bar,B_array_simu,'*r',predicted_N,predicted_B,'-b');
% xlim([1.6 2]);
% %ylim([0 20]);
% a1=xlabel('$\bar{N}$');
% %a2=ylabel('$B$');
% fig_name='./figure/N_B.eps';
% figurePostTreat


figureParameter
f1=plot(N_bar,os_fre_h,'*r',predicted_N,predicted_w,'-b');
xlim([0 5]);
%ylim([0.6 0.9]);
%ylim([0.2 0.6]);
%set(gca,'YTICK',[0.6 0.7 0.8 0.9]);
%set(gca,'YTICK',[0 0.4 0.8]);
a1=xlabel('$\bar{N}$');
a2=ylabel('$\omega$');
h1=legend("Simulation","Theory");
set(h1,'location','southwest')
fig_name='./figure/N_w.eps';
figurePostTreat
