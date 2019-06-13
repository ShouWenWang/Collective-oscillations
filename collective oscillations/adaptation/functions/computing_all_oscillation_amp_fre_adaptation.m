function [N_bar,os_fre,A_array_simu,B_array_simu]=computing_all_oscillation_amp_fre_adaptation(N_array,para_array,which_model,plot_on)

%clear all; close all;
fresh_computation=1;
if fresh_computation
        disp("Running fresh computation!!!")
%% first part, numerically determine the oscillation amplitude and frequency for the coupled system
%N_array=1:0.05:3; % range of the effective coupling used inthe paper
%N_array=[1,1.2]; % range of the effective coupling used inthe paper


%% parameters
dt=0.01;
T0=10000;
cutoff_time=500;
T_min=5; % the minimum period, this is for adaptive system
alpha1=para_array(6);
alpha2=para_array(7);

%% other initialization
B_array_simu=zeros(length(N_array),1);
A_array_simu=zeros(length(N_array),1);
rela_unc_mean_fre_h=zeros(length(N_array),1);
rela_unc_mean_amp_h=zeros(length(N_array),1);
rela_unc_mean_fre_x=zeros(length(N_array),1);
rela_unc_mean_amp_x=zeros(length(N_array),1);

os_fre_h=zeros(length(N_array),1);
os_fre_x=zeros(length(N_array),1);
os_fre_h2=zeros(length(N_array),1);
heat_array=zeros(length(N_array),1);


%% simulation
for j=1:length(N_array)
 switch which_model
     case 1
       disp("linear negative-feedback adaptation + nonlinear signal")
       wh_nonlinear=1; wx_nonlinear=0;
       [time,h,x_bar,heat_rate_h]=generating_collective_motion_feedback_adaptation(T0,dt,N_array(j),para_array,wh_nonlinear,wx_nonlinear,plot_on);
     case 2
       disp("nonlinear negative-feedback adaptation")
       wh_nonlinear=0; wx_nonlinear=1;
       [time,h,x_bar,heat_rate_h]=generating_collective_motion_feedback_adaptation(T0,dt,N_array(j),para_array,wh_nonlinear,wx_nonlinear,plot_on);
     case 3
        disp("nonlinear negative-feedback adaptation")
       [time,h,x_bar,heat_rate_h]=generating_collective_motion_feedforward_adaptation(T0,dt,N_array(j),para_array,plot_on);
 end
heat_array(j)=mean(heat_rate_h);


%[rela_unc_mean_fre_h(j),rela_unc_mean_amp_h(j),os_fre_h(j),B_array_simu(j)]=oscillation_amplitude_fre_estimation_above_zero(h,dt,T0,T_min);
[rela_unc_mean_fre_h(j),rela_unc_mean_amp_h(j),os_fre_h(j),B_array_simu(j)]=oscillation_amplitude_fre_estimation(h,dt,T0,T_min,plot_on);
%[rela_unc_mean_fre_x(j),rela_unc_mean_amp_x(j),os_fre_x(j),A_array_simu(j)]=oscillation_amplitude_fre_estimation_above_zero(x_bar,dt,T0,T_min);
[rela_unc_mean_fre_x(j),rela_unc_mean_amp_x(j),os_fre_x(j),A_array_simu(j)]=oscillation_amplitude_fre_estimation(x_bar,dt,T0,T_min,plot_on);

%%%%%%%%% an alternative way to determine the oscillation frequency
[Ch_fre,Ch_velo,fre]=Correlation_spectrum(h,dt,cutoff_time,plot_on);
%[Ch_fre,Ch_velo,fre]=Correlation_spectrum(h,dt,cutoff_time);
[amp_h,index_h]=max(Ch_fre(2:end));
os_fre_h2(j)=fre(1+index_h);

end
N_bar=N_array*alpha1*alpha2;
os_fre=os_fre_h2; %using the result determined from the correlation spectrum
signal_amp=alpha2*B_array_simu;  % renormalized amplitude 


%% save all data
% tt=today;
% aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
% output_name="./data/"+aa+"_collective_dynamics_amp_fre_"+num2str(floor(10*rand()));
output_name="./data/collective_dynamics_amp_fre";
save(output_name);




else
    
    load  ./data/collective_dynamics_amp_fre
    
end

if plot_on

    figureParameter
    f1=plot(N_bar(2:end),os_fre(2:end),'*r');%,N_bar,os_fre_x,'-*b',N_bar,os_fre_h2,'-*g');
    a1=xlabel('$\bar{N}$');
    ylim([0, 0.3]);
    xlim([0 3]);
    a2=ylabel('$\omega$');
    %h1=legend('New method:$h\;\;$','New method:$x\;\;$','Old method:$h\;$');
    %set(h1,'location','southeast')
    %ylim([0 0.3]);
    fig_name='./figure/N_w.eps';
    figurePostTreat

    figureParameter
    f1=plot(N_bar,signal_amp,'*r');%,N_bar,os_fre_x,'-*b',N_bar,os_fre_h2,'-*g');
    a1=xlabel('$\bar{N}$');
    %ylim([0, 0.6]);
    xlim([0 3]);
    a2=ylabel('$B$');
    fig_name='./figure/N_B.eps';
    figurePostTreat

end