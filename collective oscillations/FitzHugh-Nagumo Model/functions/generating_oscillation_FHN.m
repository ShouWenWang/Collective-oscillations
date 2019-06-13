function [h,aver_x,heat_rate_h]=generating_oscillation_FHN(T0,dt,epsilon,N_eff,plot_on)

%% 3 variable, feedback adaptation, spontaneous oscillation,  

%% parameters
L=floor(T0/dt);
N=1000; % number of Dyctostilum
alpha_2=1;
alpha_1=N_eff/(N*alpha_2);


station_activity=-1.513;% only when the cell activity is larger than this number, will it secretes signal!

wy=0.2; 
%epsilon=0.01; %0.5
a0=1.5;%0.8
T=0.1;


wh=1; %signal degradation rate
beta_2=0; %basal signal generation rate
beta_0=0;
%self-degradation rate
%alpha_1=1; %signal release rate
gamma_h=1;
gamma_x=1;
gamma_y=1;


%% initialization
h=zeros(L,1);
%h(floor(100/dt):end)=1;
x=zeros(L,N);
x(1,:)=randn(1);
y=zeros(L,N);
y(1,:)=randn(1);


%noise_h=randn(L,1)*sqrt(dt*2*T/gamma_h);
noise_h=zeros(L,1);
noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);
noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);

aver_x=zeros(L,1);

for j=2:L
    
    aver_x(j-1)=sum(x(j-1,:))/N;
    h(j)=h(j-1)+N*beta_2*dt+N*alpha_1*(aver_x(j-1)-station_activity)*dt-(wh+N*beta_0)*h(j-1)*dt-h(j-1)^3*dt;%+noise_h(j-1);%-h(j-1)^3*dt; %noise cannot be added, otherwise negative signal emerges, bad for log
    
    x(j,:)=x(j-1,:)+(x(j-1,:)-(x(j-1,:)).^3/3)*dt-y(j-1,:)*dt+noise_x(j-1,:)+alpha_2*h(j-1)*dt;
    y(j,:)=y(j-1,:)+wy.*(x(j-1,:)-epsilon.*y(j-1,:)+a0)*dt+noise_y(j-1,:);
end
   aver_x(L)=sum(x(end,:))/N;
N_bar=N*alpha_2*alpha_1;



%% energetics
b=1/2; % not relevant for the average heat flux

Qh=0; %Qy=0; Qx1=0; Qx2=0;
for j=2:L
    Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1));
   % Qx1=Qx1+(b*y(j)+(1-b)*y(j-1))*(x(j)-x(j-1))*wx*gamma_x;
   % Qx2=Qx2+(x(j)-x(j-1))*(b*h(j)+(1-b)*h(j-1));
   % Qy=Qy-(b*x(j)+(1-b)*x(j-1))*(y(j)-y(j-1))*wy*gamma_y;
end

heat_rate_h=alpha_1*Qh/T0;
%heat_rate_x=(Qx1+Qx2)/T0;
%heat_rate_y=Qy/T0;

%eta=heat_rate_h/(heat_rate_h+heat_rate_x+heat_rate_y);



if plot_on
%N_bar=N;
time=(1:floor(T0/dt))*dt;
if N>=2
    figureParameter
    f1=plot(time,x(:,1),'-b',time,x(:,2),'-g',time,h,'-r');
    a1=xlabel('time');
    xlim([100 300]);
    ylim([-3 5 ]);
    h1=legend('$s$','$a_1$','$a_2$');
    legend boxoff
    fig_name='./figure/cAMP-oscillation-noise_1.eps';
    figurePostTreat;
else 
    
    figureParameter
    f1=plot(time,x(:,1),'-b',time,h,'-r');
    a1=xlabel('time');
    xlim([100 300]);
    ylim([-3 5 ]);
    h1=legend('$s$','$a_2$');
    legend boxoff
    fig_name='./figure/cAMP-oscillation-noise_1.eps';
    figurePostTreat;
end

% x21=x;y21=y;
% 
% figureParameter
% f1=plot(x11(100/dt:end,1),y11(100/dt:end,1),'*r',x21(200/dt:end,1),y21(200/dt:end,1),'-k');
% a1=xlabel('$x$');
% a2=ylabel('$y$');
% %xlim([0 500]);
% %ylim([0 100 ]);
% %h1=legend('$s/N$','$a_1$','$a_2$');
% %legend boxoff
% fig_name='./figure/cAMP-oscillation-noise_2.eps';
% figurePostTreat;


% figureParameter
% f1=plot(time,h,'-r',time,aver_x,'-.g');
% a1=xlabel('time');
% xlim([0 300]);
% %ylim([-20 20]);
% h1=legend('$h$','$\bar{x}$');
% legend boxoff
% fig_name='./figure/cAMP-oscillation-noise.eps';
% figurePostTreat;

% time=(1:L)*dt;
% figureParameter
% f1=plot(time,h,'-r',time,x(:,1),'--b',time,y(:,1),'-.k');
% a1=xlabel('time');
% xlim([0 150]);
% ylim([-15 20]);
% h1=legend('$h$','$x$','$y$');
% fig_name='./figure/cAMP-oscillation-noise.eps';
% figurePostTreat;

% figure, plot(time, epsilon*g(h),'-r')
% xlabel('time');
% ylabel('a*g(h)');


% figureParameter
% f1=plot(time,x(:,1),'-r',time,x(:,2),'--b',time,x(:,3),'-.k');
% a1=xlabel('time');
% xlim([0 400]);
% %ylim([-3 6]);
% h1=legend('$x_1$','$x_2$','$x_3$');
% fig_name='./figure/cAMP-oscillation-noise.eps';
% figurePostTreat;

end



