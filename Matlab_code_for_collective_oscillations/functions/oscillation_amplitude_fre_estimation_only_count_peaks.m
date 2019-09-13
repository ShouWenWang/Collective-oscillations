function [fre_final,amp_final]=oscillation_amplitude_fre_estimation_only_count_peaks(h,dt,T0,T_min,plot_on)

%% pick the peak amplitude of the signal
%T_min  the period within which two peaks are unlikely to coexist [can be
%set to zero if not excitable system]
% it is an important parameter to filter the out the local fluctuations

% deal with the peak and trough, to eliminate the effect of vertical bias

% For this method to be accurate, T0 must be much larger than T_ini (i.e.,
% containing many periods from T_ini to T0)
%%
T_ini=3; %90, delete the transient period


L=T0/dt;
L_ini=floor(T_ini/dt);

%T_min=10;  % the minimum period
Delta=(T_min/dt)/2;

%grad_h=h(2:end)-h(1:end-1);
%curve_h=grad_h(2:end)-grad_h(1:end-1);
amp_max=zeros(10,1);
time_max=zeros(10,1);



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

%%
amp_max_mean=mean(amp_max);
amp_final=amp_max_mean-mean(h(L_ini,:));

% relative uncertainty of the mean amplitude


%%

period=time_max(2:end)-time_max(1:end-1);

fre=2*pi./period;


fre_final=mean(fre);
uncertainty_fre=std(fre);
% relative uncertainty of the mean frequency
%rela_unc_mean_fre=uncertainty_fre/(sqrt(length(fre))*fre_final);



%% criteria for oscillation

if plot_on
    if uncertainty_fre/fre_final>1.5
        %% figures
        figure,hist(amp_max,20);
        xlabel('amp');
        figure,hist(period,20);
        xlabel('period');
        error('Fluctuation of period is too large, not justified as oscillation');
    end
end





