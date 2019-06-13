
function [alpha,beta,phi,lambda,aver_X,X_index,distr_x,T]=Langevin_2d_FRR(which_model)

%which_model=1;
switch which_model
    case 1
       
%% model 1: adaptation
disp("Running nonlinear adaptation model")
%T=1
% X_0=-3;X_1=3;
% Y_0=-4;Y_1=4;
% dx=0.15;
% dy=0.25;

%%T=0.1;
% X_0=-1;X_1=1;
% Y_0=-2;Y_1=2;
% dx=0.1;
% dy=0.1;

%%T=0.01;
X_0=-0.5;X_1=0.5;
Y_0=-1;Y_1=1;
dx=0.05;
dy=0.05;

epsilon=0;

h0=0;
wx=1;
wy=1;%wy=0.5;
nonlinear=1;
gamma_x=1; gamma_y=1; T=0.01;
Ux=@(x,y) gamma_x*(wx*1/2*(x-y).^2+1/4*nonlinear*x.^4+h0*x); %
Uy=@(x,y) gamma_y*(wy*x.*y+epsilon*y.^2/2);  %F=-x;  
% to prevent instability of the algorithm

case 2
%% model 2: linear adaptation
disp("Running linear adaptation model")

X_0=-0.5;X_1=0.5;
Y_0=-0.5;Y_1=0.5;
dx=0.02;
dy=0.02;
epsilon=0;

h0=0;
wx=1;
wy=1;%wy=0.5;
gamma_x=1; gamma_y=1; T=0.01;
Ux=@(x,y) gamma_x*(wx*1/2*(x-y).^2+h0*x); %
Uy=@(x,y) gamma_y*(wy*x.*y+epsilon*y.^2/2);  %F=-x;  
% to prevent instability of the algorithm


    case 3
%% model 2:  hair bundle
disp("Running hair bundle model")

X_0=-1.8;X_1=1.8;
Y_0=-3;Y_1=3;
dx=0.1;
dy=0.25;

% 
k=2;
a=3.5;

tau=10; %10
b=4;
gamma_x=1; gamma_y=1;
T=0.1;
h=0;

 
 Ux=@(x,y) gamma_x*(1/2*k*x.^2-1/2*a*(x-y).^2+1/4*(x-y).^4-h*x);
 Uy=@(x,y) gamma_y*(-b/tau*x.*y+1/2*y.^2/tau);


    case 4
%% model 4: FNH model
disp("Running FHN model")

    which_parameter=3;
    
    switch which_parameter
        case 1 %%%%%% parameter set: 1, T=0.001
            X_0=-1.75;X_1=-1.25;
            Y_0=-0.75;Y_1=0.25;
            dx=0.025;
            dy=0.025;

            wy=0.2; 
            epsilon=0.1;
            x0=1.5;
            gamma_x=1; gamma_y=1;
            T=0.001; k=1;
            a=1;

        case 2  %%%%%% parameter set: 2, T=0.01
            X_0=-3;X_1=2.5;
            Y_0=-1.5;Y_1=2;
            dx=0.2;
            dy=0.2;

            wy=0.2; 
            epsilon=0.1;
            x0=1.5;
            gamma_x=1; gamma_y=1;
            T=0.01; k=1;
            a=1;

        case 3 %%%%%% parameter set: 3, T=0.1
            X_0=-3;X_1=3;
            Y_0=-2;Y_1=4;
            dx=0.3;
            dy=0.3;


            wy=0.2; 
            epsilon=0.1;
            x0=1.5;
            gamma_x=1; gamma_y=1;
            T=0.1; k=1;
            a=1;
    end

    Ux=@(x,y) gamma_x*(-k/2*x.^2+a*1/12*x.^4+y.*x); 
    Uy=@(x,y) gamma_y*(-wy*x.*y+1/2*epsilon*wy*y.^2-wy*x0*y);


    case 5
%% cAMP, log sensing
disp("Running cAMP log-sensing model")

% X_0=-3;X_1=3;
% Y_0=-1;Y_1=3;
% dx=0.3;
% dy=0.3;

X_0=-3;X_1=3;
Y_0=-2;Y_1=4;
dx=0.25;
dy=0.2;

wy=0.1; 
epsilon=0.5;
x0=1.2;
gamma_x=1; gamma_y=1;
T=0.1;

h0=0;
%N=1; alpha2=0.1; beta1=0.1; beta0=0.1; w2=1;
%h0=*alpha0/(J+alpha_PDE*rho);

Ux=@(x,y) gamma_x*(-1/2*x.^2+1/12*x.^4+y.*x-h0*x);
Uy=@(x,y) gamma_y*(-wy*x.*y+1/2*epsilon*wy*y.^2-wy*x0*y);


   
    case 6
%% 2-d simple model
disp("Running ? model")
wx=0.5; wy=0.3;w1=-1;
X_0=-5;X_1=5;
Y_0=-5;Y_1=5;
dx=0.2;
dy=0.2;

gamma_x=1;gamma_y=1;T=1;
Ux=@(x,y) gamma_x*(wx*1/2*x.^2-w1*x.*y);
Uy=@(x,y) gamma_y*(wy*1/2*y.^2-x.*y);

    case 7
%% signal dynamics
disp("Running ? model")
wx=1; wy=0.00001;w1=1;
X_0=-1;X_1=1;
Y_0=-1;Y_1=1;
dx=0.02;
dy=0.1;

gamma_x=1;gamma_y=1;T=0.01;
Ux=@(x,y) gamma_x*(wx*1/2*x.^2+w1*x.^4/4);
Uy=@(x,y) wy*1/2*y.^2;


    case 8
%% Van del Vol Oscillator (m=1),  treat it as 2-d overdamped langevin, gamma_x=m;
disp("Running Van del Pol model")
X_0=-5;X_1=5;
Y_0=-3;Y_1=3;
dx=0.2;
dy=0.2;
% 


lambda=1; mu=2; A=1; k=1;
gamma_x=1;gamma_y=1; T_x=1; T_y=0.01;
% m=gamma_x in our mapping

% x is the velocity, y is the position
Ux=@(x,y) (lambda-mu*(1-(y/A).^2))*x.^2/2+k*x.*y;
Uy=@(x,y) -x.*y;

T=T_x;

end

%%
if ~exist('T_x','var')
T_x=T;
T_y=T;
end
%% 
X_index=X_0:dx:X_1;
Y_index=Y_0:dy:Y_1;
%X_index=dx*(1:round(X_length/dx))-0.5*X_length;  % date point: N+1  
Nx=length(X_index);
%Y_index=dy*(1:round(Y_length/dy))-0.5*Y_length;
Ny=length(Y_index);
tot_N=Nx*Ny;


dUx=zeros(Nx-1,Ny);
dUy=zeros(Nx,Ny-1);

epsilon=0.000001; % to avoid dUx=0 or dUy=0;
for j=1:Ny
    dUx(:,j)=Ux(X_index(2:end),Y_index(j))-Ux(X_index(1:end-1),Y_index(j))+epsilon*rand(1,Nx-1); % size: Nx-1: Ny
end

for i=1:Nx
    dUy(i,:)=Uy(X_index(i),Y_index(2:end))-Uy(X_index(i),Y_index(1:end-1))+epsilon*rand(1,Ny-1); % size: Nx: Ny-1
end




b_x=dUx./(exp(dUx/T_x)-1).*1/(gamma_x*dx^2); % birth rate
d_x=-dUx./(exp(-dUx/T_x)-1).*1/(gamma_x*dx^2);% death rate

b_y=dUy./(exp(dUy/T_y)-1).*1/(gamma_y*dy^2); % birth rate
d_y=-dUy./(exp(-dUy/T_y)-1).*1/(gamma_y*dy^2);% death rate


%% transition rate matrix  (reflectory boundary condition)
    A=zeros(tot_N,tot_N);
    
for i=2:Nx-1
    for j=2:Ny-1
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-d_x(i-1,j)-b_y(i,j)-d_y(i,j-1);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp-Nx)=b_y(i,j-1);
    A(tp,tp+Nx)=d_y(i,j);
    end
end

j=1;
for i=2:Nx-1
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-d_x(i-1,j)-b_y(i,j);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp+Nx)=d_y(i,j);

end

j=Ny;
for i=2:Nx-1
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-d_x(i-1,j)-d_y(i,j-1);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp-Nx)=b_y(i,j-1);
end


i=1;    
for j=2:Ny-1
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-b_y(i,j)-d_y(i,j-1);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp-Nx)=b_y(i,j-1);
    A(tp,tp+Nx)=d_y(i,j);
end

i=Nx;    
for j=2:Ny-1
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-d_x(i-1,j)-b_y(i,j)-d_y(i,j-1);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp-Nx)=b_y(i,j-1);
    A(tp,tp+Nx)=d_y(i,j);
end

% boundary terms
i=1; j=1;
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-b_y(i,j);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp+Nx)=d_y(i,j);
    
i=1;j=Ny;
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-b_x(i,j)-d_y(i,j-1);
    A(tp,tp+1)=d_x(i,j);
    A(tp,tp-Nx)=b_y(i,j-1);

i=Nx;j=1;    
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-d_x(i-1,j)-b_y(i,j);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp+Nx)=d_y(i,j);
    
i=Nx; j=Ny;   
    tp=Nx*(j-1)+i; 
    A(tp,tp)=-d_x(i-1,j)-d_y(i,j-1);
    A(tp,tp-1)=b_x(i-1,j);
    A(tp,tp-Nx)=b_y(i,j-1);

%% eigenvalues,  distribution,  dissipation
[NormVector,orderEigValue]=orderedEigSystem(A,1);

Prob=NormVector(:,1)';

norm_Prob=zeros(Nx,Ny);
for j=1:Ny
   tp=Nx*(j-1);
   norm_Prob(1:Nx,j)=Prob(tp+1:tp+Nx); 
end

% figure, surf(X_index,Y_index,norm_Prob');
% xlabel('x');
% ylabel('y');
% zlabel('Prob');

figureParameter
f1=pcolor(X_index,Y_index,norm_Prob');
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/distribution.eps';
figurePostTreat


distr_x=zeros(1,Nx);
for i=1:Nx
    distr_x(i)=sum(norm_Prob(i,:));
end
distr_x=distr_x./sum(distr_x)*1/dx;

distr_y=zeros(1,Ny);
for i=1:Ny
    distr_y(i)=sum(norm_Prob(:,i));
end
distr_y=distr_y./sum(distr_y)*1/dy;

figure, plot(X_index,distr_x,'-r');
xlabel('x');
ylabel('Prob-x');



%% FDR  (correlation and response),  coefficients (alpha, beta, phi)

Obser_X=zeros(1,tot_N);



for j=1:Ny
    tp=Nx*(j-1)+1;
    Obser_X(tp:tp+Nx-1)=X_index;
end

aver_X=sum(Prob.*Obser_X);

alpha=zeros(1,tot_N);
for i=1:tot_N
    alpha(i)=sum(Obser_X.*NormVector(:,i)');
    %alpha(i)=sum(Obser_X.^2.*NormVector(:,i)');
end


% it is very important to have conjugate
Y=conj(inv(NormVector));
%Y=inv(NormVector);
beta=zeros(1,tot_N);
for i=1:tot_N
    beta(i)=sum(Y(i,:).*Obser_X.*Prob);
    %beta(i)=sum(Y(i,:).*Obser_X.^2.*Prob);
end


B=zeros(1,tot_N);   %perturbation vector
for i=2:Nx-1
    for j=1:Ny
    tp=Nx*(j-1)+i; 
    B(tp)=dx/(2*T)*((Prob(tp-1)*b_x(i-1,j)+Prob(tp)*d_x(i-1,j))-(Prob(tp)*b_x(i,j)+Prob(tp+1)*d_x(i,j)));
    end
end

    i=1;
    for j=1:Ny
    tp=Nx*(j-1)+i; 
    B(tp)=dx/(2*T)*(-(Prob(tp)*b_x(i,j)+Prob(tp+1)*d_x(i,j)));
    end
    
    i=Nx;
    for j=1:Ny
    tp=Nx*(j-1)+i; 
    B(tp)=dx/(2*T)*((Prob(tp-1)*b_x(i-1,j)+Prob(tp)*d_x(i-1,j)));
    end

%%

phi=zeros(1,tot_N);
for i=1:tot_N
    phi(i)=sum(B.*Y(i,:));
end


lambda=-orderEigValue(:,1)'; %eigenvalues.


%% check the sum rule
sumrule=0;
for i=2:tot_N
    sumrule=sumrule+alpha(i)*(phi(i)-beta(i)*lambda(i));
end

