%%  phase shift of cell and signal for nonlinear adaptive circuits
tau_h=0.1;
tau_a=1;tau_y=1;
Rh_fre1=@(sigma,tau_h) 1./(1-sqrt(-1)*sigma*tau_h);
gamma_x=1;
Rx_fre_epsilon=@(sigma,epsilon) 1/gamma_x*1./(1-sqrt(-1)*sigma*tau_a-1./(sqrt(-1)*sigma*tau_y-epsilon));
G=@(sigma,B)  Rh_fre1(sigma,tau_h).*Rx_fre_epsilon(sigma,epsilon);

%omega=0.00001:0.03:100;
%omega1=0.00001:0.001:100;
omega=0.001.*1.3.^(0:100);
omega1=0.001.*1.1.^(0:100);
%Rh=Rh_fre1(omega);


figureParameter
f1=semilogx(omega,angle(Rh_fre1(omega,10)),'*r',omega,angle(Rh_fre1(omega,1)),'og',omega,angle(Rh_fre1(omega,0.1)),'sk',omega1,-angle(Rx_fre_epsilon(omega1,0.1)),'-b');
a1=xlabel('$\omega$');
xlim([0.001 100]);
ylim([0,1.6]);
set(gca,'Ytick',[0 1.6]);
set(gca,'Xtick',[ 0.01  1  100 ]);
a2=ylabel('$\phi_a^*$');
% h1=legend('$\tau_s=10$','$\tau_s=1$','$\tau_s=0.1$');%,'$\epsilon=0.8$');
% %h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
% set(h1,'location','northwest')
fig_name='./figure/omega_angleR1.eps';
figurePostTreat


figureParameter
f1=semilogx(omega,angle(Rh_fre1(omega,10)),'*r',omega1,-angle(Rx_fre_epsilon(omega1,0.01)),'-k',omega1,-angle(Rx_fre_epsilon(omega1,0.1)),'-b',omega1,-angle(Rx_fre_epsilon(omega1,0.5)),'-cyan');
a1=xlabel('$\omega$');
xlim([0.001 100]);
ylim([0,1.6]);
set(gca,'Ytick',[0 1.6]);
set(gca,'Xtick',[ 0.01  1  100 ]);
a2=ylabel('$\phi_a^*$');
% h1=legend('$\phi_s$','$\epsilon=0.01$','$\epsilon=0.1$','$\epsilon=0.5$');%,'$\epsilon=0.8$');
% %h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
% set(h1,'location','northwest')
fig_name='./figure/omega_angleR2.eps';
figurePostTreat



%% numerical simulation


T0=10000;dt=0.01;


N=10;
wx_nonlinear=1; wh_nonlinear=0;
plot_on=1;

tau_s=1;
tau_a=1;
tau_y=1;
gamma_a=1;
epsilon=0.1;
alpha1=0.5;
alpha2=0.5;
para_array=[tau_s,tau_a,tau_y,gamma_a,epsilon,alpha1,alpha2];

[time,h,x_bar,heat_rate_h]=generating_collective_motion_feedback_adaptation(T0,dt,N,para_array,wh_nonlinear,wx_nonlinear,plot_on);





%% compute the numerically exact correlation and response spectrum
   Omega=0.001.*1.1.^(0:300);
   Time=0.01*(1:10000);%0.0001.*1.1.^(0:100);
   
which_model=1;
[alpha,beta,phi,lambda,aver_X,X,Prob_x,T]=Langevin_2d_FRR(which_model);
output_cell=correlation_response(alpha,beta,phi,lambda);
R_time=output_cell{1,1};
R_fre_velo=output_cell{4,1};
C_time=output_cell{5,1};
C_fre_velo=output_cell{8,1};
FRR_vio_velo_fre=output_cell{9,1};
Chi_time=output_cell{10,1};
FRR_vio_integral=output_cell{11,1};

data_cell={R_time(Time),R_fre_velo(Omega),C_time(Time),...
    C_fre_velo(Omega),FRR_vio_velo_fre(Omega),Chi_time(Time),FRR_vio_integral};

k=1;



%%%%%%%%%%%%  test the FRR in high frequency limit

figureParameter;
f1=semilogx(Omega, real(data_cell{k,4}),'-or',Omega,2*T*real(data_cell{k,2}),'*k');%fre,Cx_velo,'--b');%Omega,2.*Rv(Omega),'-.k');
%xlim([0.01,100]);
xlim([0.001,100]);
set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
%ylim([-0.1,1]);
a1=xlabel('$\omega$');
%a2=ylabel('$\tilde{C}_{\dot{x}}(\omega)-2T\tilde{R}_{\dot{x}}''(\omega)$');
h1=legend('$ \tilde{C}_{v}$','$2T\tilde{R}_{v}''$');%,'$ \tilde{C}_{v}^*$');
set(h1,'location','northwest');
fig_name='./figure/FRR_violation_cAMP_T1.eps';
figurePostTreat


figureParameter;
f1=semilogx(Omega,-imag(data_cell{k,2})./Omega,'-r',Omega,real(data_cell{k,2})./Omega,'*k',Omega, 0.5*real(data_cell{k,4})./(T*Omega),'ob');%wr_1,zeros(1,length(index)),'+g');
xlim([0.01,100]);
%ylim([-2,12]);
set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
%ylim([-1,15]);
%set(gca, 'YTick' , [-2 0 2 4 6 8 10 12] );
a1=xlabel('$\omega$');
%a2=ylabel('$\tilde{C}_{\dot{x}}(\omega)-2T\tilde{R}_{\dot{x}}''(\omega)$');
h1=legend('$\tilde{R}_{x}''$','$ \tilde{R}_{x}''''$','$ \frac{1}{2T\omega}\tilde{C}_{x}$');
set(h1,'location','northeast');
%legend boxoff
fig_name='./figure/FRR_violation_cAMP_T1.eps';
figurePostTreat

%%%%%%%%%%% response in the time domain

Chi_x=cumsum(R_time(Time)).*(Time(2)-Time(1));

time1=-10:0.1:100;
y=zeros(1,length(time1));
y(100:end)=1;

figureParameter
f1=plot(time1,y,'-r',Time,Chi_x-1.53,'-b');
xlim([-10 40]);
ylim([-2 2])
a1=xlabel('Time');
h1=legend('$s$','$\langle a_t\rangle$');
set(h1,'location','east');
legend boxoff
fig_name='./figure/excitable_system.eps';
figurePostTreat


%%%%%%%%%%%%% response amplitude and phase shift
Omega1=0.001.*1.2.^(0:100);
wh=1; alpha1=1;
Rh_fre=@(sigma) 1./(wh-sqrt(-1)*sigma);
Rx_fre=@(sigma) sqrt(-1).*R_fre_velo(sigma)./sigma;
Cx=C_fre_velo(Omega1)./Omega1.^2;



figureParameter, 
f1=semilogx(Omega1,-angle(Rx_fre(Omega1)),'-b');%Omega,angle(Rh_fre(Omega)),'-.b');
a1=xlabel('$\omega$');
xlim([0.01,100]);
set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
ylim([-pi/2,pi/2]);
a2=ylabel('$\phi_a$');
set(gca,'YTICK',[-1.6 0 1.6]);
%h1=legend('$\phi_a$','$-\phi_s$');
set(h1,'location','southwest');
%set(gca,'YTick',[-pi/2 0 pi/2]);
fig_name='./figure/angle_adaptation_signal_T1.eps';
figurePostTreat

figureParameter, 
f1=semilogx(Omega1,abs(Rx_fre(Omega1)),'-b');
a1=xlabel('$\omega$');
a2=ylabel('$|\tilde{R}_a|$');
xlim([0.01,100]);
set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
set(gca,'YTick',[0 1 2]);
fig_name='./figure/angle_adaptation_amp.eps';
figurePostTreat


%% testing an old theory about the spectrum of the joined system, that is no longer very interesting

%wh=3*10^(-4);gamma_h=1;T=0.001;epsilon=1;
wh=3.5;gamma_h=1;T=0.1;alpha_1=1; alpha_2=1;
Rx_fre=alpha_2*sqrt(-1)*data_cell{k,2}./Omega;
Ch0_fre=(2*T/gamma_h)./(wh^2+Omega.^2);
Rh_fre=alpha_1*1/gamma_h*1./(wh-sqrt(-1)*Omega);
Cx=data_cell{k,4}./Omega.^2;
Rv_real=real(data_cell{k,2});
Cv_fre=data_cell{k,4};
Ch_theory_2=(Ch0_fre+abs(Rh_fre).^2.*Cx)./(abs(1-Rh_fre.*Rx_fre)).^2;
J_fre=(Ch0_fre/(2*T)).*(2*T*Rv_real-Cv_fre)./(abs(1-Rh_fre.*Rx_fre)).^2;


figureParameter, semilogx(Omega,Ch_theory_2,'-r');
xlim([10^(-2) 10]);
%ylim([0 1]);
a1=xlabel('$\omega$');
fig_name='./figure/cAMP_Ch_theory.eps';
figurePostTreat

figureParameter
f1=semilogx(Omega,Ch_theory_2,'-r',fre,real(Ch_fre),'+b');
xlim([10^(-2) 10]);
%ylim([0 1]);
a1=xlabel('$\omega$');
h1=legend('$\tilde{C}_h^*$:theory','$\tilde{C}_h^*$:simu');
set(h1,'location','northeast');
fig_name='./figure/cAMP_Ch_theory_simulation.eps';
figurePostTreat


%figure, semilogx(Omega,Ch_theory_2,'-r');

figureParameter 
f1=semilogx(Omega,real(Rx_fre),'-r',Omega,imag(Rx_fre),'*r',Omega,real(1./Rh_fre),'-g',Omega,imag(1./Rh_fre),'+g');
a1=xlabel('$\omega$');
xlim([0.01 100]);
ylim([-1 1]);
h1=legend('$\tilde{R}_x''$','$\tilde{R}_x''''$','$(1/\tilde{R}_h)''$','$(1/\tilde{R}_h)''''$');
figurePostTreat


figureParameter 
f1=semilogx(Omega,(real(1-Rh_fre.*Rx_fre)).^2,'+k',Omega,(imag(1-Rh_fre.*Rx_fre)).^2,'*r',Omega,(real(1-Rh_fre.*Rx_fre)).^2,'-g');
a1=xlabel('$\omega$');
xlim([0.01 100]);
ylim([0.00001 10]);
h1=legend('$[(1-\tilde{R}_h\tilde{R}_x)'']^2$','$[(1-\tilde{R}_h\tilde{R}_x)'''']^2$','sum');
set(h1,'location','southeast');
figurePostTreat


%%  null cline
%cAMP

 x0=1.5;a=1;wh=1;
f=@(x,N) x-x.^3/3+N*a*x/wh;

g=@(x,epsilon) 1/epsilon*(x+x0);
x=-20:0.1:20;
figureParameter 
f1=plot(x,f(x,1),'sr',x,f(x,5),'+k',x,f(x,10),'*b',x,g(x,0.1),'-g',x,g(x,0.01),'--g');
a1=xlabel('$x$');
a2=ylabel('$y$');
xlim([-8 8]);   
ylim([-30 30]);   
%h1=legend('$N=1$','$N=15$','$N=30$');
%set(h1,'location','southeast');
fig_name='./figure/cAMP_Ch_nullcline_epsilon_001.eps';
figurePostTreat



%%  nonlinear response
%wr=0.1678;
h0=0.1*2.^(1:10);
which_model=2;
for i=1:length(h0)
[Omega,Rx_fre,Rv_fre]=nonlinear_response(h0(i),which_model);
data_nonlinear_response_dicty{i,1}=Omega;
data_nonlinear_response_dicty{i,2}=Rx_fre;
data_nonlinear_response_dicty{i,3}=Rv_fre;
end


N_fre=length(data_nonlinear_response_dicty{i,2});
Rx_real_array=zeros(N_fre,length(h0));
Rx_imag_array=zeros(N_fre,length(h0));
Omega1=data_nonlinear_response_dicty{1,1};
for i=1:length(h0)
Rx_real_array(:,i)=real(data_nonlinear_response_dicty{i,2});
Rx_imag_array(:,i)=imag(data_nonlinear_response_dicty{i,2});
end

figureParameter
loglog(h0,Rx_real_array(1,:),'-+r',h0,-Rx_imag_array(1,:),'-*b');
xlim([0.1 200])
ylim([0.001 1])
set(gca,'xtick',[0.001 0.01 0.1 1 10 100]);
a1=xlabel('$h_0$');
h1=legend('$\tilde{R}_x''$','$\tilde{R}_x''''$');
set(h1,'location','southwest');
fig_name='./figure/dicty-nonlinear_response_h0_01.eps';
figurePostTreat

figureParameter
loglog(h0,Rx_real_array(4,:),'-+r',h0,-Rx_imag_array(4,:),'-*b');
xlim([0.1 200])
ylim([0.001 1])
set(gca,'xtick',[0.001 0.01 0.1 1 10 100]);
a1=xlabel('$h_0$');
h1=legend('$\tilde{R}_x''$','$\tilde{R}_x''''$');
set(h1,'location','southwest');
fig_name='./figure/dicty-nonlinear_response_h0_16.eps';
figurePostTreat

%[Omega,Rx_fre,Rv_fre]=nonlinear_response_adaptation(h0);
% Omega=0.001.*1.2.^(0:100);
figureParameter
f1=semilogx(Omega1,Rx_real_array(:,1),'-+r',Omega1,Rx_imag_array(:,1),'-+r',...
    Omega1,Rx_real_array(:,2),'-sk',Omega1,Rx_imag_array(:,2),'-sk',...
    Omega1,Rx_real_array(:,3),'-.b',Omega1,Rx_imag_array(:,3),'-.b',...
    Omega1,Rx_real_array(:,4),'-og',Omega1,Rx_imag_array(:,4),'-og',...
    Omega1,Rx_real_array(:,5),'--r',Omega1,Rx_imag_array(:,5),'--r',...
    Omega1,Rx_real_array(:,6),'-+k',Omega1,Rx_imag_array(:,6),'-+k',...
    Omega1,Rx_real_array(:,7),'-*b',Omega1,Rx_imag_array(:,7),'-*b',...
    Omega1,Rx_real_array(:,8),'-.g',Omega1,Rx_imag_array(:,8),'-.g',...
    Omega,real(data_cell{k,2})./Omega,'-cyan',Omega,-imag(data_cell{k,2})./Omega,'-cyan');
xlim([0.01 100])
set(gca,'xtick',[0.001 0.01 0.1 1 10 100]);
a1=xlabel('$\omega$');
%a2=ylabel('$\tilde{R}_x''$');
%h1=legend('$\tilde{R}_x'''':1$','$\tilde{R}_x'':1$','$\tilde{R}_x'''':0$','$\tilde{R}_x'':0$');
%set(h1,'location','northwest');
fig_name='./figure/dicty-nonlinear_response.eps';
figurePostTreat


%% phase diagram
%%T=1
data1=[[0.2 4];[0.1 10];[0.08 13];[0.04 28];[0.02 61];[0.01 129]];
%T=0.1
data=[[0.4 1.8];[0.3 2.5];[0.2 4.5];[0.1 10];[0.05 20];[0.02 60];[0.01 129]];


figureParameter
loglog(data(:,1),data(:,2),'-*r',data(:,1),1.2+zeros(1,length(data(:,1))),'-*r');
xlim([0.01 1]);
a1=xlabel('$\epsilon$');
a2=ylabel('$\bar{N}$');
fig_name='./figure/Dicty_phase_diagram.eps';
figurePostTreat

%% phase and amplitude

wh=1;
Rh_fre=1./(wh-sqrt(-1)*Omega);
Rx_fre=sqrt(-1)*data_cell{k,2}./(Omega);
phi_h=-angle(Rh_fre);
phi_x=-angle(Rx_fre);

% figureParameter, 
% f1=semilogx(Omega,abs(Rx_fre),'ob',Omega,1./(abs(Rh_fre)),'-r',Omega,1./(2*abs(Rh_fre)),'-.k',Omega,1./(9*abs(Rh_fre)),'+g');%Omega,abs(Rh_fre)/5,'-+k',Omega,abs(Rh_fre)/11,'-*cyan');
% a1=xlabel('$\omega$');
% xlim([0.001,100]);
% set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
% ylim([0,1.5]);
% %h1=legend('$|\tilde{R}_x|$','$|\tilde{R}_h|$','$|\tilde{R}_h/2|$','$|\tilde{R}_h|/4$','$|\tilde{R}_h|/20$');
% %set(h1,'location','northeast');
% %set(gca,'YTick',[-pi/2 0 pi/2]);
% fig_name='./figure/amplitude_adaptation_signal_T1.eps';
% figurePostTreat


%w_r=[[3 0.153];[4 0.115];[5 0.096];[6 0.081];[7 0.067];[8 0.058];[9 0.048];[10 0.04];[11 0.033]];
w_r=[[1 0.2];[2 0.1246];[3 0.1];[ 4 0.08149];[5 0.07191];[6 0.05752];[7 0.04794];[8 0.038]];
figureParameter, 
f1=semilogx(Omega,-angle(Rx_fre),'-.r',Omega,angle(Rh_fre),'-b',w_r(:,2),zeros(1,length(w_r(:,2))),'*g');
a1=xlabel('$\omega$');
xlim([0.001,100]);
set(gca, 'XTick' , [0.001 0.01 0.1 1 10 100 1000] );
ylim([-pi/2,pi/2]);
set(gca,'YTICK',[-1.57 0 1.57]);
h1=legend('$\phi_x$','$-\phi_h$');
set(h1,'location','southwest');
%set(gca,'YTick',[-pi/2 0 pi/2]);
fig_name='./figure/angle_adaptation_signal_T1.eps';
figurePostTreat

%% resonance rates
% % T=1;
% resonance_data=[3 0.153 0.1534  0.1534
% 4 0.115  0.115  0.1246
% 5 0.096  0.096   0.1055
% 6  0.081  0.08149   0.08629
% 7 0.067  0.0671   0.0719
% 8  0.058  0.058   0.0599
% 9 0.048  0.048 0.0503
% 10  0.04  0.038  0.04075
% 11 0.033  0.02876  0.02876];

resonance_data=[
1  0.1678  0.19  0.2
2  0.1246    0.1246  0.1342
3  0.1    0.1  0.1055
4 0.08149   0.08149  0.08629
5  0.07191   0.07191  0.0767
6  0.05752   0.05752  0.05752
7  0.04794  0.04794  0.04314
8  0.038    0.038 0.039];

figureParameter, 
f1=plot(resonance_data(:,1),resonance_data(:,2),'-or');
a1=xlabel('$\bar{N}$');
a2=ylabel('$\omega$');
xlim([1 8]);
set(gca, 'XTick' , [1 2 3 4  5 6 7 8]); %6 8 10 12] );
ylim([0,0.2]);
fig_name='./figure/resonance_N.eps';
figurePostTreat



figureParameter, 
f1=plot(resonance_data(:,1),resonance_data(:,2),'ob',...
    resonance_data(:,1),resonance_data(:,3),'--r',...
    resonance_data(:,1),resonance_data(:,4),'+k',resonance_data(:,1),0.3./resonance_data(:,1),'-g');
a1=xlabel('$\bar{N}/\omega_2$');
a2=ylabel('$\omega$');
xlim([1 8]);
set(gca, 'XTick' , [1 2 3 4  5 6 7 8]); %6 8 10 12] );
ylim([0,0.2]);
h1=legend('$w_2=1,\alpha_1=1$','$w_2=1,\alpha_1=0.5\;\;$','$w_2=2,\alpha_1=1$');
set(h1,'location','northeast');

fig_name='./figure/resonance_N.eps';
figurePostTreat


%% eigenvalue

 x0=1.5;wh=1; epsilon=0.08;kappa=1;
 dx=0.01;

N_temp=1:0.5:15;
largest_eig=zeros(length(N_temp),1);
second_eig=zeros(length(N_temp),1);
for k=1:length(N_temp);%length(N_temp)
N=N_temp(k);

%sol_x=vpasolve(x-x.^3/3+N*a*A*x/wh-(x+x0)/epsilon==0,x);
Poly=[-1/3 0 1+N/wh-1/epsilon -x0/epsilon];
sol_x=roots(Poly);
j=1;
while abs(imag(sol_x(j)))>0.0001 && j<3
    j=j+1;
end

sol_x_real=real(sol_x(j));
g=@(x) x-x^3/3+N*x/wh;%-(x+x0)/epsilon;
kx=(g(sol_x_real+dx)- g(sol_x_real-dx))/(2*dx);
A=[[kx -1];  [kappa -kappa*epsilon]];
[NormVector,orderEigValue]=orderedEigSystem(A);

largest_eig(k)=orderEigValue(1,1);
second_eig(k)=orderEigValue(2,1);
end

figureParameter
f1=plot(N_temp,real(largest_eig),'-or',N_temp,imag(largest_eig),'-+b');
a1=xlabel('$\bar{N}$');
a2=ylabel('$\lambda_{max}$');
h1=legend('real','imag');
set(h1,'location','southwest');
fig_name='./figure/largest_eigenvalue.eps';
figurePostTreat

figureParameter
f1=plot(N_temp,real(second_eig),'-or',N_temp,imag(second_eig),'-+b');
a1=xlabel('$\bar{N}$');
a2=ylabel('$\lambda_{2}$');
ylim([-5 1]);
h1=legend('real','imag');
set(h1,'location','southwest');
fig_name='./figure/eigenvalue_2.eps';
figurePostTreat

%% nullcline

 x0=0.9; epsilon=0.1;  

 x=-2:0.01:2;
 y1=x-x.^3/3;
 y2=(x+x0)/epsilon;
 
figure,plot(x,y1,'-r',x,y2,'-.b') 
ylim([-0.68 -0.64]);
%xlim([-1.05 -0.95]);
 
figure,plot(x(1:end-1),y1(2:end)-y1(1:end-1),'-r')
ylim([-0.1 0.1]);

%% velocity field plot

 x0=1.5;a=1;wh=1; kappa=1; epsilon=0.1;  N=20; 
max_x=100;
[x,y] = meshgrid(-max_x:1:max_x,-max_x:1:max_x);


u =(N+1)*x-x.^3/3-y;
v =kappa*(x+x0-epsilon*y);
u=0.1*u./sqrt(u.^2+v.^2);
v=0.1*v./sqrt(u.^2+v.^2);

figureParameter
f1=quiver(x,y,u,v);
a1=xlabel('$x$');
a2=ylabel('$y$');
%xlim([-6 -3]);
%ylim([-50 -30]);
xlim([-3 -2]);
ylim([-30 -8]);
%set(gca,'YTICK',[-50 -30 -10 10 30 50 ]);
set(gca,'YTICK',[-15  -10 -5 0 5 10 15 20]);% 10 30 50 ]);
set(gca,'XTICK',[-6 -4 -2 0 2 4 6  ]);
fig_name='./figure/velocity_field.eps';
figurePostTreat

%%
epsilon=0.1; dt=0.01; N=1; T0=1000;
[x1,y1,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=1; T0=1000;
[x2,y2,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=1; T0=1000;
[x3,y3,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=1; T0=1000;
[x4,y4,time]=reverse_motion(T0,dt,N,epsilon);

figureParameter
f1=plot(x1(:,2),y1(:,2),'-b',x3(:,1),y3(:,1),'-k',x2(:,3),y2(:,3),'-b',x2(:,2),y2(:,2),'-r',x2(:,1),y2(:,1),'-r');%x4(:,1),y4(:,1),'-b');%x2(:,1),y2(:,1),'-r',x2(:,2),y2(:,2),'-r',x3(:,1),y3(:,1),'-b',x3(:,2),y3(:,2),'-b');%x3(:,2),y3(:,2),'-r');%x4(:,3),y4(:,3),'-r',x4(:,2),y4(:,2),'-r');%x5(:,1),y5(:,1),'-r',x5(:,2),y5(:,2),'-r');
a1=xlabel('$x$');
a2=ylabel('$y$');
xlim([-5 5]);
ylim([-5 5]);
fig_name='./figure/velocity_field.eps';
figurePostTreat


%%
epsilon=0.1; dt=0.01; N=5; T0=1000;
[x1,y1,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=5; T0=1000;
[x2,y2,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=5; T0=1000;
[x3,y3,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=5; T0=1000;
[x4,y4,time]=reverse_motion(T0,dt,N,epsilon);


figureParameter
f1=plot(x2(:,1),y2(:,1),'-k',x3(:,1),y3(:,1),'-r',x1(:,1),y1(:,1),'-k',x1(:,2),y1(:,2),'-r',x1(:,1),y1(:,1),'-r');%x4(:,1),y4(:,1),'-b');%x2(:,1),y2(:,1),'-r',x2(:,2),y2(:,2),'-r',x3(:,1),y3(:,1),'-b',x3(:,2),y3(:,2),'-b');%x3(:,2),y3(:,2),'-r');%x4(:,3),y4(:,3),'-r',x4(:,2),y4(:,2),'-r');%x5(:,1),y5(:,1),'-r',x5(:,2),y5(:,2),'-r');
a1=xlabel('$x$');
a2=ylabel('$y$');
xlim([-10 10]);
ylim([-15 15]);
fig_name='./figure/velocity_field.eps';
figurePostTreat


%% phase portrait
epsilon=0.1; dt=0.01; N=15; T0=1000;
[x1,y1,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=15; T0=1000;
[x2,y2,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=15; T0=1000;
[x3,y3,time]=reverse_motion(T0,dt,N,epsilon);

epsilon=0.1; dt=0.01; N=15; T0=1000;
[x4,y4,time]=reverse_motion(T0,dt,N,epsilon);


figureParameter
f1=plot(x1(:,1),y1(:,1),'-k',x3(:,2),y3(:,2),'-r',x3(:,1),y3(:,1),'-r');%x4(:,1),y4(:,1),'-b');%x2(:,1),y2(:,1),'-r',x2(:,2),y2(:,2),'-r',x3(:,1),y3(:,1),'-b',x3(:,2),y3(:,2),'-b');%x3(:,2),y3(:,2),'-r');%x4(:,3),y4(:,3),'-r',x4(:,2),y4(:,2),'-r');%x5(:,1),y5(:,1),'-r',x5(:,2),y5(:,2),'-r');
a1=xlabel('$x$');
a2=ylabel('$y$');
xlim([-15 15]);
set(gca,'xtick',[-15 -10 -5 0 5 10 15]);
ylim([-60 60]);
fig_name='./figure/velocity_field.eps';
figurePostTreat

x2(:,1),y2(:,1),'-b',x4(:,1),y4(:,1),'-b',
%x3(:,3),y3(:,3),'-k'





epsilon=0.1; dt=0.01; N=5; T0=1000;
[x4,y4,time]=reverse_motion(T0,dt,N,epsilon);
