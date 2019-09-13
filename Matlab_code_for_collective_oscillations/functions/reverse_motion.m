function [x,y,time]=reverse_motion(T0,dt,N,epsilon,which_model)

%% 3 variable, feedback adaptation, spontaneous oscillation,  


%dt=0.01;
L=floor(T0/dt);

switch which_model
    case 1   %%  adaptation model
        kappa=0.01; 
        u=@(x,y) (N-1)*x-x.^3+y;
        v=@(x,y) -kappa*x-kappa*epsilon*y;
    
    case 2  %%  FN model
        kappa=0.2; x0=1.5;
        u =@(x,y) (N+1)*x-x.^3/3-y;
        v =@(x,y) kappa*(x+x0-epsilon*y);
end

x=zeros(L,3);
y=zeros(L,3);


% x(1,1)=-5; %initial condition
% y(1,1)=-100;
% 
% x(1,2)=5; %initial condition
% y(1,2)=100;
% 
% x(1,3)=-5.2; %initial condition
% y(1,3)=-40;


% x(1,1)=3;
% y(1,1)=9.2;
% 
x(1,1)=-2;
y(1,1)=-9.9;
% 
% x(1,3)=-3;
% y(1,3)=-9;
% 
% x(1,1)=4;
% y(1,1)=42.9;

% 
% x(1,1)=-20;
% y(1,1)=-43.3;


for j=2:L
    x(j,:)=x(j-1,:)+u(x(j-1,:),y(j-1,:))*dt; 
    y(j,:)=y(j-1,:)+v(x(j-1,:),y(j-1,:))*dt;
end
    
time=(1:L)*dt;



% figure, plot(time,x(:,1),'-r',time,x(:,2),'-.k',time,x(:,3),'--b',time,x(:,4),'-g');
% xlabel('time');
% xlim([0 500]);
% ylabel('x_bar');

% 
% figureParameter
% f1=plot(time,h/N,'-r',time,x(:,1),'-.g',time,x(:,2),'--b');
% a1=xlabel('time');
% xlim([0 50]);
% ylim([-4 6]);
% h1=legend('$h/N$','$x_1$','$x_2$');
% fig_name='./figure/Adaptation-collective-noise.eps';
% figurePostTreat;


figureParameter
f1=plot(x(:,1),y(:,1),'-.k');%x(:,2),y(:,2),'-+r',x(:,3),y(:,3),'-og');
a1=xlabel('$x$');
a2=ylabel('$y$');
xlim([-10 10]);
%ylim([-10 15]);
fig_name='./figure/Adaptation-reverse.eps';
figurePostTreat;
