function glycolysis_reaction_network
%% outdated, not used


T0=60;dt=0.005;
L=floor(T0/dt);

time=dt*(1:L)';


n=2;


S1=zeros(L,n);
S2=zeros(L,n);
S3=zeros(L,n);
S4=zeros(L,n);%+0.2; S4(50/dt:100/dt)=0.3;  
%S4=0.2+0.1*cos(1*time);
S4ext=zeros(L,1); %S4ext(50/dt:100/dt)=0.2;  
A3=zeros(L,n);
N2=zeros(L,n);

S1(1,:)=5.8;
S2(1,:)=0.9;
S3(1,:)=0.2;
S4(1,:)=0.2;%*(1+rand(1,n));
N2(1,:)=0.1;
A3(1,:)=2.4;
S4ext(1,:)=0.1;

J0=3; 
k1=100;
k2=6;
k3=16;
k4=100;
k5=1.28;
k6=12;
k=100;
% kappa=13;
kappa=13;
q=4;
Ki=0.52;
N=1;
A=4;
phi=0.1;

f=@(x) 1./(1+(x/Ki).^q);

for i=2:L
    v1=k1*S1(i-1,:).*A3(i-1,:).*f(A3(i-1,:));
    v2=k2*S2(i-1,:).*(N-N2(i-1,:));
    v3=k3*S3(i-1,:).*(A-A3(i-1,:));
    v4=k4*S4(i-1,:).*N2(i-1,:);
    v5=k5*A3(i-1,:);
    v6=k6*S2(i-1,:).*N2(i-1,:);
    v7=k*S4ext(i-1);
    Ja=kappa*(S4(i-1,:)-S4ext(i-1));
    
    S1(i,:)=S1(i-1,:)+(J0-v1)*dt;
    S2(i,:)=S2(i-1,:)+(2*v1-v2-v6)*dt;
    S3(i,:)=S3(i-1,:)+(v2-v3)*dt;
    S4(i,:)=S4(i-1,:)+(v3-v4-Ja)*dt;
    N2(i,:)=N2(i-1,:)+(v2-v4-v6)*dt;
    A3(i,:)=A3(i-1,:)+(-2*v1+2*v3-v5)*dt;
    S4ext(i)=S4ext(i-1)+(phi*sum(Ja)/n-v7)*dt;
    
end


% figureParameter
% f1=plot(time,S1(:,1),'--r',time,S2(:,1),'-.g',time,S3(:,1),'+k',...
%     time,S4(:,1),'-r',time,S4ext,'-g',time,N2(:,1),'-b',time,A3(:,1),'ocyan');
% a1=xlabel('time');
% xlim([50 60]);
% %h1=legend('ACA$_{int}$','ACA$_{ext}$','NADH');%'$x_1$','$x_2$');%,'$x_2$','$y_2$');
% h1=legend('S1','S2','S3','S4int','S4ext','N2','A3');
% set(h1,'location','northwest');
% fig_name='./figure/ACA_NADH.eps';
% figurePostTreat;




figureParameter
f1=plot(time,S4(:,1),'-*r',time,S4(:,1),'-b',time,S4ext,'-g',time,N2(:,1),'ok');
a1=xlabel('time');
xlim([50 55]);
h1=legend('ACA$_{int}1$','ACA$_{int}2$','ACA$_{ext}$','NADH');%'$x_1$','$x_2$');%,'$x_2$','$y_2$');
%h1=legend('S4int','S4ext','N2');
set(h1,'location','northwest');
fig_name='./figure/ACA_NADH.eps';
figurePostTreat;



figureParameter
f1=plot(time,S4(:,1)-mean(S4(:,1)),'-*r',time,5*(S4ext-mean(S4ext)),'-g',time,-N2(:,1)+mean(N2(:,1)),'ok',time,S3(:,1)-mean(S3(:,1)),'-k',time,0.3*(-A3(:,1)+mean(A3(:,1))),'+g');
a1=xlabel('time');
xlim([50 55]);
h1=legend('ACA$_{int}$','ACA$_{ext}$','NAD','S3','ADP');%'$x_1$','$x_2$');%,'$x_2$','$y_2$');
%h1=legend('S4int','S4ext','N2');
set(h1,'location','northwest');
fig_name='./figure/ACA_NADH.eps';
figurePostTreat;



