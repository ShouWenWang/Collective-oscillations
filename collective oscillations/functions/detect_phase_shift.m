function [theta,x_amp,h_amp]=detect_phase_shift(fre,h,dt,x)
%% we use fre, amp, dt, T0 to define h
% x is the output from simulation, and we try to determine its phase shift
% from h. 

T0=floor(length(x)*dt);
max_period=floor(T0/(2*pi/fre)-2);
if max_period<3
    error("please increase the duration")
end
%% the intermediate period for estimating the phase shift
T1=(max_period-2)*2*pi/fre; 
T2=max_period*2*pi/fre;

d=dt*fre/pi;% step 

%%
index_1=floor(T1/dt);
index_2=floor(T2/dt);
x_used=x(index_1:index_2);
x_mean=mean(x_used);
h_used=h(index_1:index_2);
h_mean=mean(h_used);
x_amp=sqrt(sum((x_used-x_mean).^2)*d)/sqrt(2);
h_amp=sqrt(sum((h_used-h_mean).^2)*d)/sqrt(2);
CosTheta=d*sum((h_used-h_mean).*(x_used-x_mean)/(h_amp*x_amp))/2;

index_quarter_period=floor(pi/(2*fre*dt));
index_3=index_1+index_quarter_period;
index_4=index_2+index_quarter_period;
h_Cos=h(index_3:index_4);
SinTheta=d*sum((h_Cos-h_mean).*(x_used-x_mean)/(h_amp*x_amp))/2;

if SinTheta>=0
    theta=acos(CosTheta);
else
    theta=-acos(CosTheta);
end
    
%% 
