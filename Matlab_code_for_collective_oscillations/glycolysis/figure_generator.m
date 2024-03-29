
startup

clear all; close all;


%% plot all variables at a given cell density, with explicit diffusion, slow
T0=10; %default 10, maximum number of period that you want to use to compute frequency and amplitude
D=100;
dt=0.000001;
plot_fig=1;
density=0.01;
[s_in_auto,s_ex_auto,time]=minimum_glycolysis(dt,T0,density,D,plot_fig);


%% plot all variables at a given cell density, eliminate diffusion 
T0=10; %default 20, maximum number of period that you want to use to compute frequency and amplitude
dt=0.00001;
plot_fig=1;
density=0.121;
[s,a,time]=minimum_glycolysis_no_Diffusion(dt,T0,density,plot_fig);


%% plot the oscillation amplitude and frequency as the cell density is changed, for Fig. 
%dt=0.01;
T0=10; %default 20, maximum number of period that you want to use to compute frequency and amplitude
T_min=0.1; % originally 0.1
%step_size=0.0001;
plot_on=1;
density=cat(2,[0.001:0.02:1,1]); 
%density=0.001:0.02:0.4; 
with_explicit_diffusion=1;
D=10;
dt=0.0001; % for D=100, you need to use dt=0.000002 to avoid blowing-up solution
Time_ini=2; % initial time when you start to estimate the period

fre_final_a=zeros(1,length(density));
amp_final_a=zeros(1,length(density));

fre_final_s=zeros(1,length(density));
amp_final_s=zeros(1,length(density));

ave_level_s=zeros(1,length(density));
ave_level_s_ex=zeros(1,length(density));
ave_level_a=zeros(1,length(density));


plot_fig=1;
L_mid=floor(2/dt);

parfor j=1:length(density)
    
if with_explicit_diffusion
    [s_in_auto,s_ex_auto,a_auto,time]=minimum_glycolysis(dt,T0,density(j),D,plot_fig);
    ave_level_s_ex(j)=mean(s_ex_auto(L_mid:end));
else
    [s_in_auto,a_auto,time]=minimum_glycolysis_no_Diffusion(dt,T0,density(j),plot_fig);
end

ave_level_a(j)=mean(a_auto(L_mid:end));
ave_level_s(j)=mean(s_in_auto(L_mid:end));


[rela_unc_mean_fre,rela_unc_mean_amp,fre_final_s(j),amp_final_s(j)]=oscillation_amplitude_fre_estimation(s_in_auto,dt,T0,T_min,plot_fig,Time_ini);
[rela_unc_mean_fre,rela_unc_mean_amp,fre_final_a(j),amp_final_a(j)]=oscillation_amplitude_fre_estimation(a_auto,dt,T0,T_min,plot_fig,Time_ini);


%[fre_final_s(j),amp_final_s(j)]=oscillation_amplitude_fre_estimation_only_count_peaks(s_in_auto,dt,T0,T_min,plot_fig);
%[fre_final_a(j),amp_final_a(j)]=oscillation_amplitude_fre_estimation_only_count_peaks(a_auto,dt,T0,T_min,plot_fig);
disp(num2str(density(j)));
end


%save data/glycolysis_all_D100 fre_final_s amp_final_s fre_final_a amp_final_a ave_level_s ave_level_a ave_level_s_ex
%load data/glycolysis_all_D01
figureParameter
left_color = [124 0 0]/255;
right_color = [0 128 0 ]/255;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
f2=plot(density,ave_level_s,'-');
ylim([0,2])
%ylabel("ACE_in")
%ylabel('')

yyaxis right
f1=plot(density,ave_level_s_ex,'-.');
fig_name='./figure/test.eps';
ylim([0,2])
xlabel("Cell density")
%ylabel("ACE_ex")
figurePostTreat



figureParameter
f1=plot(density,amp_final_s,'-b');
ylim([0 0.1])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
xlabel("Cell density")
%ylabel('Amplitude');
fig_name="./figure/density_amp_glycolysis_s.eps";
figurePostTreat

figureParameter
f1=plot(density,fre_final_s,'-b');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
xlabel("Cell density")
ylabel('Oscillation frequency-s');
fig_name="./figure/density_fre_glycolysis_s.eps";
figurePostTreat

figureParameter
f1=plot(density,ave_level_s,'-b',density,ave_level_s_ex,'-r');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
ylim([0 2])
xlabel("Cell density")
ylabel('Average signal');
fig_name="./figure/density_concentration_s.eps";
figurePostTreat


figureParameter
f1=plot(density,amp_final_a,'-b');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
xlabel("Cell density")
ylabel('Oscillation amplitude-a');
fig_name="./figure/density_amp_glycolysis_a.eps";
figurePostTreat

figureParameter
f1=plot(density,fre_final_a,'-b');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
xlabel("Cell density")
ylabel('Oscillation frequency-a');
fig_name="./figure/density_fre_glycolysis_a.eps";
figurePostTreat

figureParameter
f1=plot(density,ave_level_a,'-b');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
xlabel("Cell density")
ylabel('Average activity');
fig_name="./figure/density_concentration_a.eps";
figurePostTreat

