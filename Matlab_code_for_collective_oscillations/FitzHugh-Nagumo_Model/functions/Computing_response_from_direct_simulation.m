function [Rx_fre,Rv_fre,time,h,x_bar,Energy_outflow,efficiency,Internal_power]=Computing_response_from_direct_simulation(Hz,h0,which_model,plot_on)
%% purpose
% computing the numerical response, and also the energy flow rate, under a
% periodic perturbation


%% master parameters for all sub models  
T0=10000;
dt=0.01;
L=floor(T0/dt);
N=100; %cell number

time=(0:L-1)*dt;
%Hz=0.1; % close to the optimum frequency
%h0=1; Hz=0.2;
h=h0*cos(2*pi*Hz*time);%+1*cos(2*pi*0.6*time);%%h0*cos(2*pi*Hz*time);
%2*cos(2*pi*0.2*time);%+
% which_model=2; % which model to use

%% figures control

plot_multi_curve=0; 
plot_energy=0;  
plot_dynamics=0;

%% initialization
x=zeros(L,N);
y=zeros(L,N);
x_bar=zeros(L,1);
x(1,:)=0;


switch which_model
    case 1  
        disp("model 1: (nonlinear) sensory adaptation");
        nonlinear_coeff=1;  % whether the system is linear or not
        w0=1;w1=1; epsilon=0.1; 
        gamma_x=1;gamma_y=1;
        T=0.01;

        noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);
        noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);

        for j=2:L
            x_bar(j-1)=1/N*sum(x(j-1,:));
            %h(j)=h(j-1)+basal*dt-wh(j-1)*h(j-1)*dt+alpha_2*N*x_bar(j-1)/gamma_h*dt+noise_h(j-1)-wh_nonlinear*h(j-1)^3*dt;
            x(j,:)=x(j-1,:)-w0*(x(j-1,:)-y(j-1,:))*dt+h(j-1)/gamma_x*dt+noise_x(j-1,:)-x(j-1,:).^3*dt;
            y(j,:)=y(j-1,:)-w1.*x(j-1,:)*dt-w1.*epsilon.*y(j-1,:)*dt+noise_y(j-1,:);
        end

        x_bar(L)=sum(x(L,:))/N;

        %% energetics
    %     b=1/2; % not relevant for the average heat flux
    % 
    %     Qy=0; Qx1=0; Qx2=0;
    %     for j=2:L
    %         %Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1)); % positive if energy flows from x to h
    %         Qx1=Qx1+(b*y(j,:)+(1-b)*y(j-1,:)).*(x(j,:)-x(j-1,:))*w0*gamma_x;
    %         Qx2=Qx2+(x(j,:)-x(j-1,:))*(b*h(j)+(1-b)*h(j-1));  % positive if energy flows from h to x
    %         Qy=Qy-(b*x(j,:)+(1-b)*x(j-1,:)).*(y(j,:)-y(j-1,:))*w1*gamma_y;
    %     end
    % 
    % 
    %     Energy_outflow=-Qx2/T0;
    %     Internal_power=(Qx1+Qy)/T0;
    %     efficiency=Energy_outflow./Internal_power;

          Energy_outflow=0;
          Internal_power=0;
          efficiency=0;
    

    case 2  
        disp("model 2: FitzHugh-Nagumo model");
        %% parameter for the Van del Pol oscillator
        %gamma_x=1; k=1; T=1; w0=k/gamma_x; c=1; wc=c/gamma_x;
        %w1=1; x0=0; epsilon=0; gamma_y=1; a=1;

        % parameters for the FN model in our paper
        gamma_x=1; k=1; T=0.1; w0=k/gamma_x; c=1; wc=c/gamma_x;
        w1=0.2; x0=1.5; epsilon=0.1; gamma_y=1; a=1;

        noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);
        noise_y=randn(L,N)*sqrt(dt*2*T/gamma_y);

        for j=2:L
            x_bar(j-1)=1/N*sum(x(j-1,:));
            %h(j)=h(j-1)+basal*dt-wh(j-1)*h(j-1)*dt+alpha_2*N*x_bar(j-1)/gamma_h*dt+noise_h(j-1)-wh_nonlinear*h(j-1)^3*dt;
            x(j,:)=x(j-1,:)+w0*x(j-1,:)*dt - wc*y(j-1,:)*dt/gamma_x+h(j-1)/gamma_x*dt+noise_x(j-1,:)-a*x(j-1,:).^3*dt/(3*gamma_x);
            y(j,:)=y(j-1,:)+w1.*(x(j-1,:)+x0)*dt-w1*epsilon.*y(j-1,:)*dt+noise_y(j-1,:);
        end

        x_bar(L)=sum(x(L,:))/N;


    %     % energetics
    %     b=1/2; % not relevant for the average heat flux
    % 
    %     Qy=0; Qx1=0; Qx2=0;
    %     for j=2:L
    %         %Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1)); % positive if energy flows from x to h
    %         Qx1=Qx1-(b*y(j,:)+(1-b)*y(j-1,:)).*(x(j,:)-x(j-1,:));
    %         Qx2=Qx2+(x(j,:)-x(j-1,:))*(b*h(j)+(1-b)*h(j-1));  % positive if energy flows from h to x
    %         Qy=Qy+(b*x(j,:)+(1-b)*x(j-1,:)).*(y(j,:)-y(j-1,:))*w1*gamma_y;
    %     end
    % 
    %     Energy_outflow=-Qx2/T0;
    %     Internal_power=(Qx1+Qy)/T0;
    %     efficiency=Energy_outflow./Internal_power;

          Energy_outflow=0;
          Internal_power=0;
          efficiency=0;
      
    case 3     
            disp("model 3: Van del Pol oscillator");
            T=0.1; % temperature
            mu=5; lambda=1; % friction
            m=1; k=1;
            wv1=lambda/m; wv2=mu/m; wv3=k/m;

            gamma_x=1;
            v=zeros(L,N);
            noise_v=randn(L,N)*sqrt(dt*2*T*lambda)/m;
            noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);


            for j=2:L
                x_bar(j-1)=1/N*sum(x(j-1,:));
                %h(j)=h(j-1)+basal*dt-wh(j-1)*h(j-1)*dt+alpha_2*N*x_bar(j-1)/gamma_h*dt+noise_h(j-1)-wh_nonlinear*h(j-1)^3*dt;
                v(j,:)=v(j-1,:)-(wv1-wv2*(1-x(j-1,:).^2)).*v(j-1,:)*dt - wv3*x(j-1,:)*dt+h(j-1)/m*dt+noise_v(j-1,:);
                % it is important that noise cannot be added here,  otherwise the
                % system is unstable
                x(j,:)=x(j-1,:)+v(j-1,:)*dt;%+noise_x(j-1,:);
            end

            x_bar(L)=sum(x(L,:))/N;


            % energetics
            b=1/2; % not relevant for the average heat flux

            Qx1=0; Qx2=0;
            for j=2:L
                %Qh=Qh+(b*x(j,:)+(1-b)*x(j-1,:))*(h(j)-h(j-1)); % positive if energy flows from x to h
                Qx1=Qx1+(b*mu*(1-x(j-1,:).^2).*v(j-1,:)+(1-b)*mu*(1-x(j,:).^2).*v(j,:)).*v(j-1,:)*dt;
                Qx2=Qx2+(x(j,:)-x(j-1,:))*(b*h(j)+(1-b)*h(j-1));  % positive if energy flows from h to x
            end

            Energy_outflow=-Qx2/T0;
            Internal_power=Qx1/T0;
            efficiency=Energy_outflow./Internal_power;

            y=v;


    case 4     
           disp("model 4: nonlinear signal dynamics, Dicty model")
           c0=0.01; c1=1; b0=0.001;  M=1; alpha2=10^4;
           h=h.*(h>0);
           for j=2:L  
              x_bar(j-1)=sum(x(j-1,:))/N; 
              x(j,:)=x(j-1,:)+dt*(c0*M-(c1+b0*M)*x(j-1,:)+M*alpha2*h(j)); 
           end
           x_bar(L)=sum(x(L,:))/N;

            Energy_outflow=0;
            Internal_power=0;
            efficiency=0;

        

    case 5 
           disp("model 5: nonlinear circuit dynamics, Dicty model");
           x0=1.2; 
           kappa=0.1; epsilon=0.5; K0=10^(-3); T=0.1;alpha1=0.18;
           gamma_x=1;

           % process of the signal
           %h=h0+h0*cos(2*pi*Hz*time)+0.001;
           %f=@(s) alpha1*log(1+s/K0);
           %h=f(h);

           noise_x=randn(L,N)*sqrt(dt*2*T/gamma_x);

           for j=2:L  
              x_bar(j-1)=sum(x(j-1,:))/N; 
              x(j,:)=x(j-1,:)+dt*(x(j-1,:)-x(j-1,:).^3/3-y(j-1,:))+noise_x(j,:)+h(j)*dt; 
              y(j,:)=y(j-1,:)+dt*kappa*(x(j-1,:)-epsilon*y(j-1,:)+x0);
           end
           x_bar(L)=sum(x(L,:))/N;

            Energy_outflow=0;
            Internal_power=0;
            efficiency=0;
    
end
%% response spectrum 
[Rx_fre,Rv_fre]=calculate_response_spectrum(x_bar,h,Hz,dt,plot_on);
[Cx_fre,Cv_fre,fre]=Correlation_spectrum(x,dt,10000,plot_on);
%% figures

if plot_dynamics==1

figureParameter
f1=plot(x(10/dt:50/dt,1),y(10/dt:50/dt,1),'-b');%time,x(:,1),'--g',time,x(:,2),'-.b');
a1=xlabel('$x$');
a2=ylabel('$y$');
%xlim([0 60]);
%ylim([-1.2 1.2]);
%h1=legend('$\omega_2$');
fig_name='./figure/phase-portrait.eps';
figurePostTreat;


figureParameter
f1=plot(time,h,'-g',time,x(:,1),'-b');%time,x(:,1),'--g',time,x(:,2),'-.b');
a1=xlabel('Time');
%a2=ylabel('$x$');
xlim([0 200]);
box(gca,'off');
%ylim([-1.2 1.2]);
%h1=legend('$\omega_2$');
fig_name='./figure/Adaptation-collective-wh.eps';
figurePostTreat;




end


%% energetics


if plot_energy==1

%%

[Prob,x]=hist(Energy_outflow,20);
Prob=Prob./sum(Prob)*1/(x(2)-x(1));

figureParameter
f1=plot(x,Prob,'-or');
a1=xlabel('Output Power');
a2=ylabel('Distribution');
%xlim([-0.5 1]);
%box(gca,'off');
%ylim([-1.2 1.2]);
%h1=legend('$\omega_2$');
fig_name='./figure/Adaptation-power-distribution.eps';
figurePostTreat;

%%
[Prob,x]=hist(Internal_power,20);
Prob=Prob./sum(Prob)*1/(x(2)-x(1));

figureParameter
f1=plot(x,Prob,'-or');
a1=xlabel('Internal Power');
a2=ylabel('Distribution');
%xlim([-0.5 1]);
%box(gca,'off');
%ylim([-1.2 1.2]);
%h1=legend('$\omega_2$');
fig_name='./figure/Adaptation-power-distribution.eps';
figurePostTreat;


%%
[Prob,x]=hist(efficiency,20);
Prob=Prob./sum(Prob)*1/(x(2)-x(1));

figureParameter
f1=plot(x,Prob,'-or');
a1=xlabel('Efficiency');
a2=ylabel('Distribution');
%xlim([-0.5 1]);
%box(gca,'off');
%ylim([-1.2 1.2]);
%h1=legend('$\omega_2$');
fig_name='./figure/Adaptation-efficiency-distribution.eps';
figurePostTreat;

end


%% plot, for multi-curves
 
if plot_multi_curve==1

[Prob_10,x_10]=hist(eta_10,20);
Prob_10=Prob_10./sum(Prob_10)*1/(x_10(2)-x_10(1));

[Prob_100,x_100]=hist(eta_100,20);
Prob_100=Prob_100./sum(Prob_100)*1/(x_100(2)-x_100(1));

[Prob_1000,x_1000]=hist(eta_1000,20);
Prob_1000=Prob_1000./sum(Prob_1000)*1/(x_1000(2)-x_1000(1));


figureParameter
f1=plot(x_10,Prob_10,'-*r',x_100,Prob_100,'-ob',x_1000,Prob_1000,'-+g');
a1=xlabel('Input power');
a2=ylabel('Distribution');
xlim([-0.2 0.4]);
%box(gca,'off');
%ylim([0 0.2]); 
set(gca, 'XTICK',[-0.2 0 0.2 0.4]);
h1=legend('$L=10^1$','$L=10^2$','$L=10^3$');
set(h1,'location','northwest');
legend boxoff
fig_name='./figure/Adaptation-efficiency-distribution.eps';
figurePostTreat;


%% variance
L_vector=[10 20 50 100 200 500 1000];
std_vector=[std(Energy_outflow_10) std(Energy_outflow_20)  std(Energy_outflow_50)  std(Energy_outflow_100)  std(Energy_outflow_200)  std(Energy_outflow_500)  std(Energy_outflow_1000)];


figureParameter
f1=loglog(L_vector,std_vector,'-or');
a1=xlabel('Time: $L$');
a2=ylabel('Power std');
%xlim([-0.5 1]);
%box(gca,'off');
%ylim([0 0.2]); 
%set(gca, 'YTICK',[0 0.1 0.2]);
%h1=legend('$L=10^1$','$L=10^2$','$L=10^3$');
legend boxoff
fig_name='./figure/Adaptation-power-std.eps';
figurePostTreat;



end


