function [time,h,x_bar,heat_rate_h]=generating_collective_motion_feedforward_adaptation(T0,dt,N,para_array,wh_nonlinear,wx_nonlinear,plot_on)
%%% simulating the collective dynamics of coupled adaptive circuits that
%%% implement negative feedback

%% parameters: normal collective oscillation

tau_h=para_array(1); wh0=1/tau_h;
tau_a=para_array(2); wx=1/tau_a;
tau_y=para_array(3); wy=1/tau_y;
gamma_x=para_array(4); 
epsilon=para_array(5);
alpha_10=para_array(6);
alpha_20=para_array(7);

L=floor(T0/dt);
%wh_nonlinear=0; % 1 for a cubic nonlinear term in h dynamics
%wx_nonlinear=1; % 1 for a cubic nonlinear term in x dynamics
basal=0; % basal signal generation rate
Temp=0.01; % temperature

gamma_h=1;
gamma_y=1;

%% condition for high-Q oscillation and entrainment 
% basal=0.5;
% Temp=0.001;
% wy=0*randn(1,N)+100;  N=100; T0=100;dt=0.001;epsilon=0.01; alpha_1=1; alpha_2=0.2;
% N=100;
% wh=zeros(L,1)+5;
% period=0.33; L_p=period/dt;
% for k=1:10
%  temp_k=54/dt+(2*k+1)*L_p;   
% wh(temp_k+1:temp_k+L_p)=0;
% end




%% heterogeneous simulation
% wy_vector=[0.25^2 0.5^2 0.75^2 1^2 1.25^2 1.5^2 1.75^2 2^2];% 2.25^2 2.5^2 2.75^2 3^2 3.25^2 3.5^2 3.75^2 4^2];
% M=length(wy_vector);
% wy(1:M)=wy_vector;
% for k=1:39
% wy(k*M+(1:M))=wy_vector;
% end
% 
% N=length(wy);
%%

%wy=0.4*(1:N);
%wy(1)=0.25^2;wy(2)=1;wy(3)=3^2;wy(4)=5^2;%;wy(5)=40;%wy(6)=40;
%alpha_1=sqrt(1); 


%% initialization
alpha_1=ones(L,1)*alpha_10;
alpha_2=ones(L,1)*alpha_20;
wh=zeros(L,1)+wh0;


h=zeros(L,1); %z=zeros(L,1);z(floor(39.5/dt))=6;
time=(1:L)*dt;
x=zeros(L,N);
y=zeros(L,N);
x_bar=zeros(L,1);

y(1,:)=rand(1,N); %initial condition
x(1,:)=rand(1,N);

%noise_h=randn(L,1)*sqrt(dt*2*Temp/gamma_h);
noise_h=zeros(L,1);
noise_x=randn(L,N)*sqrt(dt*2*Temp/gamma_x);
noise_y=randn(L,N)*sqrt(dt*2*Temp/gamma_y);


%% simulation
for j=2:L
    x_bar(j-1)=1/N*sum(x(j-1,:));
    h(j)=h(j-1)+basal*dt-wh(j-1)*h(j-1)*dt+alpha_1*N*x_bar(j-1)/gamma_h*dt+noise_h(j-1)-wh_nonlinear*h(j-1)^3*dt;
    x(j,:)=x(j-1,:)-wx*(x(j-1,:)-y(j-1,:)-alpha_2*h(j-1))*dt+noise_x(j-1,:)-0*x(j-1,:).^2*dt-wx_nonlinear*x(j-1,:).^3*dt;
    y(j,:)=y(j-1,:)-wy.*(y(j-1,:)+alpha_2*h(j-1))*dt+noise_y(j-1,:);
end
    

x_bar(L)=sum(x(L,:))/N;



%% energetics
b=1/2; % not relevant for the average heat flux

Qh=0;Qy=0; Qx1=0; Qx2=0;
for j=2:L
    Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1));
    Qx1=Qx1+(b*y(j,:)+(1-b)*y(j-1,:))*(x(j,:)-x(j-1,:))'*wx*gamma_x;
    Qx2=Qx2+(x(j,:)-x(j-1,:))*(b*h(j)+(1-b)*h(j-1));
    Qy=Qy-(b*x(j,:)+(1-b)*x(j-1,:))*(y(j,:)-y(j-1,:))'*wy*gamma_y;
end

heat_rate_h=alpha_1*mean(Qh,2)/T0; % average heat flow for individual cells
% heat_rate_x=(Qx1+Qx2)/T0;
% heat_rate_y=Qy/T0;
% 
% eta=heat_rate_h/(heat_rate_h+heat_rate_x+heat_rate_y);


%% plot figures

if plot_on==1


sample=1:20:length(time);
N_bar=N*alpha_10^2;
figureParameter
f1=plot(time,h/N_bar,'-k',time,x(:,1),'-r',time,sum(y,2)/N_bar,'-cyan');
%f1=plot(time(sample),x(sample,1),'-r',time(sample),x(sample,2),'-cyan',time(sample),x(sample,3),'-b',time(sample),x(sample,4),'-g');
a1=xlabel('time');
xlim([0 150]);
%ylim([-6 10]);
%set(gca,'ytick',[-2  0 2 4]);
h1=legend('$h/\bar{N}$','$a_1$','$y/\bar{N}$');
%legend boxoff
fig_name='./figure/Adaptation-collective-noise.eps';
figurePostTreat;


end


