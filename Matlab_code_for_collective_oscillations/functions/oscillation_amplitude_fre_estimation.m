function [rela_unc_mean_fre,rela_unc_mean_amp,fre_final,amp_final]=oscillation_amplitude_fre_estimation(h,dt,Tot_time,Time_window,plot_on,Time_0)

%% pick the peak amplitude of the signal
%Time_window  the period within which two peaks are unlikely to coexist [can be
%set to zero if not excitable system]
% it is an important parameter to filter the out the local fluctuations

% deal with the peak and trough, to eliminate the effect of vertical bias

% For this method to be accurate, Tot_time must be much larger than Time_0 (i.e.,
% containing many periods from Time_0 to Tot_time)
%%
if ~exist('Time_0')
    Time_0=90; %delete the transient period
end
L=Tot_time/dt;
L_ini=floor(Time_0/dt);

%Time_window=10;  % the minimum period
Delta=(Time_window/dt)/2;

%grad_h=h(2:end)-h(1:end-1);
%curve_h=grad_h(2:end)-grad_h(1:end-1);
amp_max=zeros(10,1);
amp_min=zeros(10,1);
time_max=zeros(10,1);
time_min=zeros(10,1);


%% pick the peak amplitude and its time
k=1;

for j=L_ini:L-3-Delta
    if h(j)>h(j-1) && h(j)>h(j+1)
      if  max(h(j-Delta:j+Delta))==h(j) % confirm that h(j) is the local maximum
        amp_max(k)=h(j);
        time_max(k)=j*dt;
        k=k+1;
      end
    end
end

%% pick the trough amplitude and its time
k=1;
for j=L_ini:L-3-Delta
    if h(j)<h(j-1) && h(j)<h(j+1)
      if  min(h(j-Delta:j+Delta))==h(j) % confirm that h(j) is the local minimum
        amp_min(k)=h(j);
        time_min(k)=j*dt;
        k=k+1;
      end
    end
end

%%
amp_max_mean=mean(amp_max);
amp_min_mean=mean(amp_min);
shift=(amp_max_mean+amp_min_mean)/2;
amp_corr=amp_max-shift;

amp_final=mean(amp_corr);
uncertainty_amp=std(amp_corr);
% relative uncertainty of the mean amplitude
rela_unc_mean_amp=uncertainty_amp/(sqrt(length(amp_corr))*amp_final);

%%

period_max=time_max(2:end)-time_max(1:end-1);
period_min=time_min(2:end)-time_min(1:end-1);
period=cat(1,period_max,period_min);

fre=2*pi./period;


fre_final=mean(fre);
uncertainty_fre=std(fre);
% relative uncertainty of the mean frequency
rela_unc_mean_fre=uncertainty_fre/(sqrt(length(fre))*fre_final);



%% criteria for oscillation

if plot_on
    if uncertainty_fre/fre_final>1.5
        %% figures
        figure,hist(amp_max,20);
        xlabel('amp');
        figure,hist(period_max,20);
        xlabel('period');
        error('Fluctuation of period is too large, not justified as oscillation');
    end
end





