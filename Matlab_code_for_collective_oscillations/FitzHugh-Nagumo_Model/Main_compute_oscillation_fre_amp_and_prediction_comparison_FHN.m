startup;

%% it is possible that this algorithm only works for alpha1=1; 

running_fresh_computation=1;
if running_fresh_computation
    
    disp("running_fresh_computation!!!")
    plot_on=0;
    
%% oscillation amplitude and frequency from directly simulating the collective dynamics
N_array=0.8:0.05:3; % range of the effective coupling used in the paper
%N_array=[0.8 0.85 1.2]; % test
[N_bar,os_fre_h,signal_amp]=computing_all_oscillation_amp_fre_FHN(N_array,plot_on);


%% predicted oscillation amplitude and frequency from numerically computed response spectrum
%%%%%%%% deciding the frequency for perturbation
% for the inferrence to be precise, we only need well-resolved spectrum
% around the frequency where oscillations actually happens [this can be
% obtained by collective simulations]
    Hz=0.03:0.002:0.05; % used in the paper
    %Hz=[0.03 0.04 0.05]; % for test,must contain 3 elements
 
  
%%%%%% deciding the amplitude of the perturbation
    h0=0.2:0.2:3; %used in the paper
    %h0=0.01:0.1:0.2;
    %h0=[0.2  0.4]; % for test
[fre_os_esti,Neff_os_esti]=predicting_amp_fre_from_computed_response_spectrum(Hz,h0,plot_on);


%% save_all_results

% tt=today;
% aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/amp_fre_predicted_comparison"+num2str(floor(10*rand()));
save(output_name);

else
    load ./data/190625_amp_fre_predicted_comparison_paper
   

end


%%% plot
figureParameter
f1=plot(N_bar(2:end),os_fre_h(2:end),'*r',Neff_os_esti(1),fre_os_esti(1),'ob');
a1=xlabel('$\bar{N}$');
ylabel('Frequency');
% xlim([0 2]);
% ylim([0 0.3]);
h=legend('Simulation','Prediction');
set(gca,'YTICK',[0 0.1 0.2 0.3]);
fig_name='./figure/N_w.eps';
figurePostTreat

figureParameter
%f1=plot(Neff_os_esti,h0,'-.b',N_bar(4:end),alpha1*B_array_simu(4:end),'*r');
f1=plot(N_bar(2:end),signal_amp(2:end),'*r',Neff_os_esti(1),h0(1),'ob');
a1=xlabel('$\bar{N}$');
ylabel('Amplitude');
h=legend('Simulation','Prediction');
%h1=legend('');
% xlim([0 2]);
%ylim([0 2]);
fig_name='./figure/N_B.eps';
figurePostTreat

