%% estimating the response specturm of dicty 
%%%% these are the data reproduced from the paper [Sgro et al. Mol. Sys. Biol. 2015], Fig. 2, at 1nM perturbation
Time=4/515*[0,31, 68, 94, 125, 160, 190, 223, 256, 287, 320, 352, 385, 419, 450, 482, 514,546,577];
FRET=[0,5, 46, 79, 104,  115, 120, 124, 116, 97, 76, 54, 35, 25, 15, 6,3,1,0];
Delta_t=0.25;
Time_q=0:Delta_t:4.5;
%%

%% short time response
figure
FRET_q = interp1(Time,FRET,Time_q,'spline');
plot(Time,FRET,'o',Time_q,FRET_q,':.');
title('Spline Interpolation');

%% add additional temporal information
additional_Time=4.5+Delta_t:Delta_t:100;
additional_FRET=zeros(1,length(additional_Time));
FRET_final=cat(2,FRET_q,additional_FRET);
Time_final=cat(2,Time_q,additional_Time);

figureParameter
semilogx(Time_final,FRET_final,'-+r');
a1=xlabel('Time');
xlim([0.1 10]);
fig_name='./figure/a_exp.eps';
figurePostTreat

%% compute the temporal response function by differentiating the above FRET signal
diff_FRET=zeros(1,length(FRET_final)-1);
diff_Time=zeros(1,length(FRET_final)-1);

for j=1:length(FRET_final)-1
    diff_FRET(j)=(FRET_final(j+1)-FRET_final(j))/(Time_final(j+1)-Time_final(j));
    diff_Time(j)=(Time_final(j+1)+Time_final(j))/2;
end

figureParameter
semilogx(diff_Time,diff_FRET,'-+r')
xlim([0.1 10]);
a1=xlabel('Time');
fig_name='./figure/Ra_exp.eps';
figurePostTreat


%% compute the response spectrum by Fourier transformation
Y=fft(diff_FRET);
L=length(diff_FRET);
figure,plot(1:length(Y),abs(Y));

P2 = Delta_t*abs(Y/L);
P1 = P2(1:L/2+1);
%P1(2:end-1) = 2*P1(2:end-1);

Fs=1/Delta_t;
f = Fs*(0:(L/2))/L;
Omega=2*pi*f;

semilogx(Omega,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
%xlabel('f (Hz)')
xlabel('Omega');
ylabel('|P1(f)|')


Ang2 = angle(Y/L);
Ang1 = Ang2(1:L/2+1);
%Ang1(2:end-1) = 2*P1(2:end-1);



%% plot the phase shift within a small frequency window
figureParameter
f1=plot(Omega,Ang1,'*-r');
%title('Single-Sided Phase shift of X(t)')
%xlabel('f (Hz)')
xlim([0.2,1])
ylim([0,1.57])
set(gca, 'YTick' , [0 1.57]);
xlabel('Omega');
ylabel('angle')
fig_name='./figure/angle_Ra_fre.eps';
figurePostTreat

%% test for effect of signal timescale
tau_s=0.35;
phi_Rs=@(w) -atan(w*tau_s);

figureParameter
f1=plot(Omega,Ang1,'*-r',Omega,-phi_Rs(Omega),'-+b');
%title('Single-Sided Phase shift of X(t)')
%xlabel('f (Hz)')
xlim([0.2,1])
ylim([0,1.57])
set(gca, 'YTick' , [0 1.57]);
xlabel('Omega');
ylabel('angle')
fig_name='./figure/angle_Ra_fre.eps';
figurePostTreat




%% Below are another example of phase shift from adaptation, just to get some intuition
% 
% t1=1; t2=2;
% integral_Response=@(t) -exp(-t/t1)+exp(-t/t2);
% Delta_t=0.1;
% time=0:Delta_t:1000;
% %figure,semilogx(time,R(time),'-r');
% 
% 
% % alias_effect=@(x) (2*pi*x).*sin(2*pi*x)./(2-2*cos(2*pi*x)); % x=fre/fre_sampling, range: 0: pi
% % %dt_simu=0.05*0.01; % this is the update time in simulation
% % fre_sp=2*pi/Delta_t; %sampling omega
% % 
% 
% 
% int_R=integral_Response(time);
% diff_R=zeros(1,length(time)-1);
% diff_t=zeros(1,length(time)-1);
% for j=1:length(time)-1
%     diff_R(j)=(int_R(j+1)-int_R(j))/(time(j+1)-time(j));
%     diff_t(j)=(time(j+1)+time(j))/2;
% end
%  
% figure,plot(diff_t,diff_R,'*r')
% xlim([0,10])
% 
% L=length(diff_R);
% n = 2^nextpow2(L);
% Y=fft(diff_R,n);
% 
% 
% 
% P2 = Delta_t*Y;
% Ra_fre = P2(1:n/2+1);
% %P1(2:end-1) = P1(2:end-1);
% 
% Fs=1/Delta_t;
% f = Fs*(0:(n/2))/n;
% Omega=2*pi*f;
% 
% figure
% semilogx(Omega,abs(Ra_fre)) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% %xlabel('f (Hz)')
% xlabel('Omega');
% ylabel('|Ra|')
% 
% 
% 
% figure,semilogx(Omega,angle(Ra_fre));
% title('Single-Sided Phase shift of X(t)')
% %xlim([0,100])
% %xlabel('f (Hz)')
% xlabel('Omega');
% ylabel('Angle')

