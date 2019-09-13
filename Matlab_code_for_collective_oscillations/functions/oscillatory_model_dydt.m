function dydt=oscillatory_model_dydt(t,y)

%% 3 variable, feedback adaptation, spontaneous oscillation,  
w11=0.5;w12=1;w22=0.5;w32=1;w33=0;
dydt=[-w11*y(1)+y(2)-y(1)^3; w12*y(1)-w22*y(2)+w22*y(3); -w32*y(2)-w33*y(3)];

% 

%% two variables, spontaneous oscillation, cilium model
%y(1): Delta
%y(2): fm
% 
% dydt=zeros(2,1);
% 
% xi=1; k=1; tau=2; alpha1=3; alpha3=1;
% dydt(1)=(y(2)-k*y(1))/xi;
% dydt(2)=-(y(2)-alpha1*dydt(1)+alpha3*dydt(1)^3)/tau;

%% two variables, spontaneous oscillation
% w1=0.1; w2=0.1;
% dydt=[-w1*y(1)+y(2)-y(1)^3;-w2*y(2)+y(1)];

%% 3 variables,  feedforward adaptation,  oscillation
% w11=1;w12=1;w22=0.1;w33=0.1;
% dydt=[-w11*y(1)+y(2); w12*y(1)-w22*y(2)-w33*y(3)-y(2)^3; y(1)-w33*y(3)];




