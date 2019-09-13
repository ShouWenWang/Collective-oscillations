  
function [Cxy_fre,fre]=Cross_correlation_spectrum(x,y,dt,cutoff_time)
%%
% obtain correlation spectrum by truncating the long time noise for
% Cx(time)

% cutoff_time,  choose to be around 10*tau_s (tau_s:  the slow timescale of the system)

%%
L2=length(x);

Cxy_time_0=ifft(fft(x,L2).*conj(fft(y,L2)))/L2; % zero padding is irrelevant since the data is essentially uncorrelated for L

Cxy_time=Cxy_time_0(1:L2/2+1);
time=dt*(1:L2/2+1);
figure, plot(time,real(Cxy_time),'-r',time,imag(Cxy_time),'-b');
xlim([0 500]);
xlabel('time');
legend('real','imag');
ylabel('Cxy');

%% Frequency spectrum from Cxy_time: truncation of long time  noise
% this truncation introduce slow modulation power spectrum in high
% frequency domain. 


if round(cutoff_time/dt)<length(Cxy_time)
    L3=round(cutoff_time/dt);
else
    L3=length(Cxy_time);
end
newL=2^nextpow2(L3);


% newTime=dt*(1:2*newL);
% newCx=zeros(2*newL,1);
% 
% cutoff=1/dt;

newCxy(1:newL)=Cxy_time_0(1:newL);
%alpha=1; %obtained by fitting to log Cx(t)
%newCx(cutoff+1:newL)=Cxy_time_0(cutoff+1).*exp(-alpha*dt*(0:newL-cutoff-1));

%newCx(1:newL)=Cxy_time_0(1:newL);
newCxy(2*newL:-1:newL+1)=Cxy_time_0(L2:-1:L2-newL+1);



%% final result
fre=2*pi*(0:newL-1)'/(2*newL*dt);
Cxy_fre=fft(newCxy)*dt;
Cxy_fre=Cxy_fre(1:newL);


%Cv_fre=fre.^2.*Cxy_fre;

