function [h,aver_x_new,output,dt]=generating_oscillation_cAMP_N(N)
%% FHN 


%% 3 variable, feedback adaptation, spontaneous oscillation,  
T0=600;
dt=0.005;
%N=4;  % number of Dyctostilum
L=floor(T0/dt);


kappa=0.1;
epsilon=0.5;
x0=1.2;
%rho=0.0001;
J=1; alpha_PDE=0.001;
w2=J+alpha_PDE*N; %signal degradation rate
c0=0.01; %basal signal generation rate
beta_2=c0*N;
%self-degradation rate

alpha_2=10000; %signal release rate
alpha_1=0.18; %0.058;


T=0.01;
%sigma=0.15;
%noise_h=randn(L,1)*sqrt(dt*2*T/gamma_h);
%noise_x=randn(L,N)*sqrt(dt)*sigma;
noise_x=randn(L,N)*sqrt(dt*2*T);
%noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);


h=zeros(L,1);
x=zeros(L,N);
y=zeros(L,N);



K0=10^(-3); 
g=@(s) log(1+s/K0);%./(1+s/K1));% signal to free energy
%f=@(x) x-x.^3/3;



%% initialization
y(1,:)=rand(1,N); %-0.45; %initial condition
x(1,:)=rand(1,N); %-1.4;

%g=@(s) s;



x_shift=0;


for j=2:L
    
    aver_x=sum((x(j-1,:)>0).*x(j-1,:))/N;
    h(j)=h(j-1)+beta_2*dt+N*alpha_2*aver_x*dt-w2*h(j-1)*dt; %noise_h(j-1);noise cannot be added, otherwise negative signal emerges, bad for log
    
    x(j,:)=x(j-1,:)+((x(j-1,:)-x_shift)-(x(j-1,:)-x_shift).^3/3)*dt-y(j-1,:)*dt+noise_x(j-1,:)+alpha_1*log(1+h(j-1)/K0)*dt;
    y(j,:)=y(j-1,:)+kappa.*(x(j-1,:)-x_shift-epsilon.*y(j-1,:))*dt+kappa.*x0.*dt; %+noise_y(j-1,:);
end

aver_x_new=zeros(1,L);
for j=1:L
   aver_x_new(j)=sum(x(j,:))/N; 
end

    
 time=(1:L)*dt;
 
%  figureParameter
%  f1=semilogy(time,h,'-r');
%  a1=xlabel('time');% a2=ylabel('$h$');
%  xlim([200 600]);
%  ylim([0.005  10^8]);
%   set(gca,'Ytick',[0.01  1  100  10000 10^6 10^8]);
%  %h1=legend('$h$');
%  fig_name='./figure/cAMP-oscillation-h.eps';
%  figurePostTreat;

 
% figureParameter
% f1=plot(time,x(:,1),'--b',time,y(:,1),'-.r',time,aver_x_new,'-g');
% a1=xlabel('time');
% xlim([200 600]);
% ylim([-3 9]);
% set(gca,'ytick',[-3 0 3 6 9]);
% h1=legend('$x$','$y$','$\bar{x}$');
% fig_name='./figure/cAMP-oscillation-x-y.eps';
% figurePostTreat;


% figure, plot(time,aver_x_new,'-g');
%  xlabel('time'); ylabel('x-mean');

output=alpha_1*g(h);

figureParameter
f1=plot(time, alpha_1*g(h),'-r',time,x(:,1),'-b',time,x(:,1),'-.g');
a1=xlabel('Time');
xlim([100 600]);
ylim([-3 6]);
set(gca,'ytick',[-3  0 3 6 ]);
%a2=ylabel('$\alpha_1 \ln (1+h/K_0)$');
fig_name='./figure/cAMP-oscillation-noise.eps';
figurePostTreat;

% figureParameter,
% f1=plot(time,x(:,1),'--b',time,x(:,2),'-.r',time,x(:,3),'-g',time,x(:,4),'-cyan');
% a1=xlabel('time');
% xlim([200 600]);
% %ylim([-3 6]);
% h1=legend('$x_1$','$x_2$','$x_3$','$x_4$');
% fig_name='./figure/cAMP-oscillation-x-sync.eps';
% figurePostTreat;


%%
% data=[
% 3   0.04794;
% 4    0.0623;
% 5    0.06711;
% 8      0.08868;
% 10  0.09108;
% 20  0.09827;
% 50  0.1007;
% 100  0.1007;
% 200  0.1031;
% 220 0.1055;
% 250 0.115;
% 300  0.1414;
% 400  0.1414;
% 500 0.1534];
% 
% figureParameter
% f1=plot(data(:,1),data(:,2),'-.+k');
% xlim([1 500]);
% a1=xlabel('Cell density: $N$');
% fig_name='./figure/cAMP-oscillation-x-sync.eps';
% figurePostTreat;

%a2=ylabel('Oscillation');
