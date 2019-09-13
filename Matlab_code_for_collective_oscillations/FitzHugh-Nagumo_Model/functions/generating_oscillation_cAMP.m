function [h,aver_x_new]=generating_oscillation_cAMP(rho)

%% 3 variable, feedback adaptation, spontaneous oscillation,  
T0=500;
dt=0.005;
N=100;  % number of Dyctostilum
L=floor(T0/dt);


epsilon=0.1;
r=0.5;
c0=1.2;
%rho=0.0001;
J=1; alpha_PDE=1000;
D=J+alpha_PDE*rho; %signal degradation rate
h0=800; %basal signal generation rate
%self-degradation rate

S=10^6; %signal release rate
a=0.058;



sigma=0.15;
%noise_h=randn(L,1)*sqrt(dt*2*T/gamma_h);
noise_x=randn(L,N)*sqrt(dt)*sigma;
%noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);


h=zeros(L,1);
x=zeros(L,N);
y=zeros(L,N);



K0=10^(-5); 
g=@(s) log(1+s/K0);%./(1+s/K1));% signal to free energy
%f=@(x) x-x.^3/3;



%% initialization
y(1,:)=rand(1,N); %-0.45; %initial condition
x(1,:)=rand(1,N); %-1.4;

%g=@(s) s;





for j=2:L
    
    aver_x=sum((x(j-1,:)>0).*x(j-1,:))/N;
    h(j)=h(j-1)+rho*h0*dt+rho*S*aver_x*dt-D*h(j-1)*dt; %noise_h(j-1);noise cannot be added, otherwise negative signal emerges, bad for log
    
    x(j,:)=x(j-1,:)+(x(j-1,:)-(x(j-1,:)).^3/3)*dt-y(j-1,:)*dt+noise_x(j-1,:)+a*log(1+h(j-1)/K0)*dt;
    y(j,:)=y(j-1,:)+epsilon.*(x(j-1,:)-r.*y(j-1,:))*dt+epsilon.*c0.*dt; %+noise_y(j-1,:);
end

aver_x_new=zeros(1,L);
for j=1:L
   aver_x_new(j)=sum(x(j,:))/N; 
end

    
 time=(1:L)*dt;
 
 figure, plot(time,h,'-r');
 xlabel('time'); ylabel('h');
 
 
 
figure,
plot(time,x(:,1),'--b',time,y(:,1),'-.r',time,aver_x_new,'-g');
xlabel('time');
xlim([0 500]);
%ylim([-3 6]);
legend('x','y','x-mean');
%fig_name='./figure/cAMP-oscillation-noise.eps';
%figurePostTreat;


% figure, plot(time,aver_x_new,'-g');
%  xlabel('time'); ylabel('x-mean');

figure, plot(time, a*g(h),'-r')
xlabel('time');
xlim([0 500]);
ylabel('a*g(h)');


% figureParameter
% f1=plot(time,x(:,1),'-r',time,x(:,2),'--b',time,x(:,3),'-.k');
% a1=xlabel('time');
% xlim([200 400]);
% ylim([-3 6]);
% h1=legend('$x_1$','$x_2$','$x_3$');
% fig_name='./figure/cAMP-oscillation-noise.eps';
% figurePostTreat;
% %% energetics
% b=1/2; % not relevant for the average heat flux
% 
% Qh=0;Qy=0; Qx1=0; Qx2=0;
% for j=2:L
%     Qh=Qh+(b*x(j)+(1-b)*x(j-1))*(h(j)-h(j-1));
%     Qx1=Qx1+(b*y(j)+(1-b)*y(j-1))*(x(j)-x(j-1))*wx*gamma_x;
%     Qx2=Qx2+(x(j)-x(j-1))*(b*h(j)+(1-b)*h(j-1));
%     Qy=Qy-(b*x(j)+(1-b)*x(j-1))*(y(j)-y(j-1))*wy*gamma_y;
% end
% 
% heat_rate_h=Qh/T0;
% heat_rate_x=(Qx1+Qx2)/T0;
% heat_rate_y=Qy/T0;
% 
% eta=heat_rate_h/(heat_rate_h+heat_rate_x+heat_rate_y);

