function [N_bar,os_fre,signal_amp]=computing_all_oscillation_amp_fre_FHN(Neff_array,plot_on)

%clear all; close all;
startup;

fresh_computation=1;
if fresh_computation
        disp("Running fresh computation!!!")
%% first part, numerically determine the oscillation amplitude and frequency for the coupled system
%Neff_array=1:0.05:3; % range of the effective coupling used inthe paper
%Neff_array=[1,1.2]; % range of the effective coupling used inthe paper


dt=0.01;
T0=1000;
%T0=1000;
% alpha2=1;  % does not work for other alpha2 for the current code
epsilon=0.1;

% other initialization
B_array_simu=zeros(length(Neff_array),1);
A_array_simu=zeros(length(Neff_array),1);
rela_unc_mean_fre_h=zeros(length(Neff_array),1);
rela_unc_mean_amp_h=zeros(length(Neff_array),1);
rela_unc_mean_fre_x=zeros(length(Neff_array),1);
rela_unc_mean_amp_x=zeros(length(Neff_array),1);

os_fre_h=zeros(length(Neff_array),1);
os_fre_x=zeros(length(Neff_array),1);
os_fre_h2=zeros(length(Neff_array),1);
heat_array=zeros(length(Neff_array),1);
cutoff_time=500;

%T_min=5; % the minimum period, this is for adaptive system
T_min=10;  % this is for FN model

for j=1:length(Neff_array)
%[time,h,x_bar,heat_rate_h]=generating_collective_motion(T0,dt,Neff_array(j),epsilon,alpha1);
%[time,h,x_bar,heat_rate_h]=generating_collective_motion_feedforward(T0,dt,Neff_array(j),alpha1);

%% 19-6-20
[h,x_bar,heat_rate_h]=generating_oscillation_FHN(T0,dt,epsilon,Neff_array(j),plot_on);




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
N_bar=Neff_array;
%os_fre=os_fre_h2; % accurate for small oscillation amplitude
os_fre=os_fre_h; % accurate for large oscillation amplitude
signal_amp=B_array_simu; 


%% save all data
% tt=today;
% aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/collective_dynamics_amp_fre_"+num2str(floor(10*rand()));
save(output_name);




else
    
    load  ./data/190406_collective_dynamics_amp_fre_1
    
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