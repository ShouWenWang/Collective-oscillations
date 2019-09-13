function [time,x,h_ext,h_int,y]=generating_collec_motion_rand_cell_signal_in_out
%% collective oscillation with spatial effect

%% parameters
% circuit
epsilon=0.1;
wx=1;
wy=1;
alpha_1=1; 

% signal
basal_int=0;
alpha_2=10; % this is related to cell density
wh_int=1; % internal signal degradation rate
wh_ext=1; % external signal degradation rate
D_0=5; % diffusion across membrane
D_1=10; % diffusion in the space
wh_nonlinear=1;

% simulation size
N=100;  % number of spatial lattices
T0=60; % simulation time
%dt=0.001;
dt=0.001;
L=floor(T0/dt);

% noise strength
T=0.01; % temperature, 
gamma_h=1;
gamma_x=1;
gamma_y=1;

% cell density control
cell_location=zeros(1,N);
density=0.7;
for l=1:density*N
    index=floor(100*rand)+1;
   cell_location(index)=1; % overlap of two identical index is possible
end

%% Simulation
h_int=zeros(L,N); %z=zeros(L,1);z(floor(39.5/dt))=6;
h_ext=zeros(L,N); %z=zeros(L,1);z(floor(39.5/dt))=6;
time=(1:L)*dt;
%h=cos(0.5*time);
x=zeros(L,N);
y=zeros(L,N);


y(1,:)=0; %initial condition
x(1,:)=0;

%noise_h=randn(L,1)*sqrt(dt*2*T/gamma_h);
noise_h=zeros(L,N)*sqrt(dt*2*T/gamma_h);
noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);
noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);



for j=2:L
    %% without conside diffusion
    x(j,:)=x(j-1,:)-wx*(x(j-1,:)-y(j-1,:))*dt+alpha_1*h_ext(j-1,:)/gamma_x*dt+noise_x(j-1,:)-0*x(j-1,:).^3*dt;
    y(j,:)=y(j-1,:)-wy.*x(j-1,:)*dt-wy.*epsilon.*y(j-1,:)*dt+noise_y(j-1,:);
    h_int(j,:)=h_int(j-1,:)+D_0*(h_ext(j-1,:)-h_int(j-1,:))*dt+basal_int*dt-wh_int*h_int(j-1,:)*dt+alpha_2*x(j-1,:).*cell_location*dt+noise_h(j-1,:)-wh_nonlinear*h_int(j-1,:).^3*dt;
    h_ext(j,:)=h_ext(j-1,:)-D_0*(h_ext(j-1,:)-h_int(j-1,:))*dt-wh_ext*h_ext(j-1,:)*dt;
    
    %% conside diffusion now
    for k=2:N-1
    h_ext(j,k)=h_ext(j,k)+D_1*(h_ext(j-1,k+1)+h_ext(j-1,k-1)-2*h_ext(j-1,k))*dt;
    end
    
    k=1;
    h_ext(j,k)=h_ext(j,k)+D_1*(h_ext(j-1,k+1)+h_ext(j-1,N)-2*h_ext(j-1,k))*dt;

    k=N;
    h_ext(j,k)=h_ext(j,k)+D_1*(h_ext(j-1,1)+h_ext(j-1,k-1)-2*h_ext(j-1,k))*dt;
   
    
end
    

% figureParameter
% f1=pcolor(time,1:N,h');
% a1=xlabel('Time');
% a2=ylabel('Space');
% fig_name='./figure/signal_spatial_evolution.eps';
% figurePostTreat

% 
figureParameter
f1=plot(time,h_ext(:,1)+20,'-g',time,h_ext(:,10)+16,'-.b',time,h_ext(:,20)+12,'-+r',time,h_ext(:,30)+8,'-*k',time,h_ext(:,40)+4,'-cyan',time,h_ext(:,50),'-+b',...
    time,h_ext(:,60)-4,'-g',time,h_ext(:,70)-8,'-.b',time,h_ext(:,80)-12,'-+r',time,h_ext(:,90)-16,'-*k',time,h_ext(:,100)-20,'-cyan',time,h_ext(:,50),'-+b');%time,x(:,1),'--g',time,x(:,2),'-.b');
a1=xlabel('Time');
ylim([-25 25]);
%a2=ylabel('All cells');
%xlim([45 78]);
%box(gca,'off');
xlim([40 60]);
%h1=legend('$\omega_2$');
fig_name='./figure/Adaptation-collective-h-ext.eps';
figurePostTreat;


figureParameter
f1=plot(time,h_ext(:,60),'-g',time,h_int(:,60),'-r');%time,x(:,1),'--g',time,x(:,2),'-.b');
a1=xlabel('Time');
%a2=ylabel('one cell');
%xlim([45 78]);
%box(gca,'off');
ylim([-2 2]);
xlim([40 60]);
h1=legend('$s_{ex}$','$s_{in}$');
fig_name='./figure/Adaptation-collective-h-in-out.eps';
figurePostTreat;

% figureParameter
% f1=plot(time,h(:,1)+30,'-g',time,h(:,3)+26,'-.b',time,h(:,5)+22,'-+r',time,h(:,7)+18,'-*k',time,h(:,9)+14,'-cyan',time,h(:,11)+10,'-+b',...
%     time,h(:,13)+6,'-.b',time,h(:,15)+2,'-+r',time,h(:,17)-2,'-*k',time,h(:,19)-6,'-cyan',time,h(:,21)-10,'-+b',...
%     time,h(:,23)-14,'-.b',time,h(:,25)-18,'-+r',time,h(:,27)-22,'-*k',time,h(:,29)-26,'-cyan',time,h(:,31)-30,'-+b');%time,x(:,1),'--g',time,x(:,2),'-.b');
% a1=xlabel('Time');
% a2=ylabel('Part cells');
% %xlim([45 78]);
% %box(gca,'off');
% xlim([10 30]);
% %h1=legend('$\omega_2$');
% fig_name='./figure/Adaptation-collective-wh2.eps';
% figurePostTreat;


figureParameter
f1=plot(1:N,h_ext(20/dt,:),'-r',1:N,h_ext(21/dt,:),'-b',1:N,h_ext(22/dt,:),'-g');
a1=xlabel('$x$');
ylim([-2 2]);
a2=ylabel('$s_{ex}$');
h1=legend('$t=20$','$t=21$','$t=22$');
set(h1,'location','southeast');
fig_name='./figure/Adaptation-space-signal.eps';
figurePostTreat


%% energetics
% b=1/2; % not relevant for the average heat flux
% 
% Qh=0;Qy=0; Qx1=0; Qx2=0;
% for j=2:L
%     Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1));
%     Qx1=Qx1+(b*y(j,:)+(1-b)*y(j-1,:))*(x(j,:)-x(j-1,:))'*wx*gamma_x;
%     Qx2=Qx2+(x(j,:)-x(j-1,:))*(b*h(j)+(1-b)*h(j-1));
%     Qy=Qy-(b*x(j,:)+(1-b)*x(j-1,:))*(y(j,:)-y(j-1,:))'*wy*gamma_y;
% end
% 
% heat_rate_h=Qh/T0;
% heat_rate_x=(Qx1+Qx2)/T0;
% heat_rate_y=Qy/T0;
% 
% eta=heat_rate_h/(heat_rate_h+heat_rate_x+heat_rate_y);

