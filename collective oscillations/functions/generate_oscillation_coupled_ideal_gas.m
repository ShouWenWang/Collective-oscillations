function [time,h,a,y]=generate_oscillation_coupled_ideal_gas


%% 
Cv=3/2;
 T=1; phi=1;V=1; N0=50; a0=6; mu0=1;
mu_gas=@(N) T*(Cv-log(V*T^Cv./(N*phi))); % chemical potential
A_gas=@(N) N.*mu_gas(N); % Helmhoze free energy, with free parameters, T, V, N
receptor=@(N,a) (N0-N-2*N0/4)*(a+a0)*mu0+100./(N0-N); % a, activity of the receptor
force=@(N) (N-N0+2*N0/4)*mu0;
N_array=0.001:0.01:N0;
%N_array=0.001:0.1:N0;

a_array=-5:0.01:5;  % array of activity
N_star=zeros(1,length(a_array));
for k=1:length(a_array)
    data=receptor(N_array,a_array(k))+A_gas(N_array);
    N_star(k)=sum(gradient(data)<0)*0.01;
end

 figureParameter
 f1=plot(a_array,N_star,'-r');
 a1=xlabel('$a$');
 a2=ylabel('$N_f^*$');
 fig_name='./figure/ideal_gas_a_N.eps';
 figurePostTreat;
 
  figureParameter
 f1=plot(a_array,force(N_star),'-r');
 a1=xlabel('$a$');
 ylim([-30 30]);
 a2=ylabel('$f_{ext}$');
 fig_name='./figure/ideal_gas_a_force.eps';
 figurePostTreat;
 
figure,plot(a_array,N_star,'-r')


figure,plot(a_array,force(N_star),'-b')



%figure, plot(N_array,mu_gas(N_array));

 figureParameter
 f1=plot(N_array,A_gas(N_array),'-r',N_array,receptor(N_array,0),'-.b',N_array,A_gas(N_array)+receptor(N_array,0),'--g');
 ylim([-100 300]);
 a1=xlabel('$N_f$');
 xlim([0 50]);
 h1=legend('$\mathcal{F}_{box}$','$\mathcal{F}_{rec}$','$\mathcal{F}$');
 set(h1,'location','Northwest');
  fig_name='./figure/ideal_gas_free_energy.eps';
 figurePostTreat;
 % 
 figure, plot(N_array,A_gas(N_array)+receptor(N_array,-1),'-b');
 ylim([20 50]);
 
%% 
T0=1000;dt=0.001;
L=floor(T0/dt);
T_rep=1;
wh=1;
wa=1;
wy=1;
epsilon=0.01;
%wyy=0.1;

alpha_1=2;

gamma_a=1;
gamma_y=1;


N_particle=zeros(L,1);
h=zeros(L,1);
a=zeros(L,1);
y=zeros(L,1);

y(1)=0; %initial condition
a(1)=0;


noise_a=randn(L,1)*sqrt(dt*2*T_rep/gamma_a);
noise_y=randn(L,1)*sqrt(dt*2*T_rep/gamma_y);



for j=2:L
    N_particle(j-1)=N_star(floor((a(j-1)+5)/0.01)+1);
    h(j-1)=force(N_particle(j-1));

    a(j)=a(j-1)-wa*(a(j-1)-y(j-1))*dt+alpha_1*h(j-1)*dt+noise_a(j-1)-(a(j-1))^3*dt;
    y(j)=y(j-1)-wy*a(j-1)*dt-wy*epsilon*y(j-1)*dt+noise_y(j-1);
   
end
    
    N_particle(L)=N_star(floor((a(L)+5)/0.01)+1);
    h(L)=force(N_particle(L));
    
%% Plots
 time=(1:L)*dt;
 figureParameter
 f1=plot(time,h);
 xlim([0 200]);
 a1=xlabel('Time');
 a2=ylabel('$h$');
 fig_name='./figure/ideal_gas_h.eps';
 figurePostTreat;
 
  figureParameter
 f1=plot(time,N_particle);
 xlim([0 200]);
 a1=xlabel('Time');
 a2=ylabel('$s$');
  fig_name='./figure/ideal_gas_N_particle.eps';
 figurePostTreat;
 
   figureParameter
 f1=plot(time,a);
 xlim([0 200]);
 a1=xlabel('Time');
 a2=ylabel('$a$');
   fig_name='./figure/ideal_gas_a.eps';
 figurePostTreat;
 
    figureParameter
 f1=plot(time,y);
 xlim([0 200]);
 a1=xlabel('Time');
 a2=ylabel('$y$');
   fig_name='./figure/ideal_gas_y.eps';
 figurePostTreat;
 
 
% figureParameter
% f1=plot(time,a,'--b',time,y,'-.k');
% a1=xlabel('time');
% %xlim([20 500]);
% %ylim([-10 15]);
% h1=legend('$a$','$y$');
% set(h1,'location','southwest');
% fig_name='./figure/Adaptation-oscillation-noise.eps';
% figurePostTreat;


%%  chemical potential