function [rela_unc_mean_fre,rela_unc_mean_amp,fre_final,amp_final]=oscillation_amplitude_fre_estimation_above_zero(h,dt,T0,T_min)

%%
% only look at the peak that is above zero. 
% this approach is particularly useful for FH model, the excitable system

% assumes that Amp_max=-Amp_min, i.e.,  symmetric oscillation around
% zero, no vertical shift


% fluctuaton of the individual sample,  and the fluctuation of the mean of
% N samples are not the same. typically,  the latter can be reduced by a
% factor of sqrt{N}



%%
T_ini=30; % delete the transient period
L=T0/dt;
L_ini=floor(T_ini/dt);

%T_min=10;  % the minimum period
Delta=floor((T_min/dt)/2);

%grad_h=h(2:end)-h(1:end-1);
%curve_h=grad_h(2:end)-grad_h(1:end-1);
amp_max=zeros(10,1);
amp_min=zeros(10,1);
time_max=zeros(10,1);
time_min=zeros(10,1);


%% pick the peak amplitude and its time
k=1;

for j=L_ini:L-3-Delta
    if h(j)>h(j-1) && h(j)>h(j+1) && h(j)>0
      if  max(h(j-Delta:j+Delta))==h(j) % confirm that h(j) is the local maximum
        amp_max(k)=h(j);
        time_max(k)=j*dt;
        k=k+1;
      end
    end
end

%% amplitude

sample_N=length(amp_max);

amp_max_mean=mean(amp_max);
amp_final=amp_max_mean;

uncertainty_amp=std(amp_max);
% relative uncertainty of the mean amplitude
rela_unc_mean_amp=uncertainty_amp/(sqrt(sample_N)*amp_final);

%% frequency

period_max=time_max(2:end)-time_max(1:end-1);

fre_max=2*pi./period_max;
fre_final=mean(fre_max);

uncertainty_fre=std(fre_max);
% relative uncertainty of the mean frequency
rela_unc_mean_fre=uncertainty_fre/(sqrt(sample_N)*fre_final);



%% criteria for oscillation
judge=0;
if uncertainty_fre/fre_final>1 && judge==1
    %% figures
    figure,hist(amp_max,20);
    xlabel('amp');
    figure,hist(period_max,20);
    xlabel('period');
    error('Fluctuation of period is too large, not justified as oscillation');
end




