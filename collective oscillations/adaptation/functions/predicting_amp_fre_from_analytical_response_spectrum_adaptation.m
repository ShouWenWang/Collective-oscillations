function [predicted_N,predicted_A,predicted_B,predicted_w]=predicting_amp_fre_from_analytical_response_spectrum_adaptation(which_model,para_array,Pert_array,plot_on)

%which_model=1;
% para_array=[tau_h,tau_a,tau_y,gamma_a,epsilon,alpha1,alpha2];

%% parameters for simulation
tau_h=para_array(1);
tau_a=para_array(2);
tau_y=para_array(3);
gamma_a=para_array(4);gamma_x=gamma_a;
epsilon=para_array(5);
alpha1=para_array(6);
alpha2=para_array(7);

%% parameters for estimating the intersection point for the phase-compensation relation
w_ini=0.001; %initial frequency
dw=0.001; % step width

%%

switch which_model

    case 1
    disp("nonlinear response theory:  linear (negative-feedback) adaptive circuit + nonlinear signal")


    Rh_fre1=@(sigma,B) 1./(1+3/4*B^2-sqrt(-1)*sigma*tau_h);%.*exp(sqrt(-1)*tau*sigma);
%     Rx_fre=@(sigma) sqrt(-1).*R_fre_velo(sigma)./sigma;
    wy=1/tau_y;
    wx=1/tau_a;
    Rx_fre_epsilon=@(sigma) 1/gamma_x*1./(wx*wy./(epsilon*wy-sqrt(-1)*sigma)-sqrt(-1)*sigma+wx);
    G=@(sigma,B)  Rh_fre1(sigma,B).*Rx_fre_epsilon(sigma);

    omega=0.00001:0.001:100;
    %omega=0.01.*1.1.^(0:100);
    Rx=Rx_fre_epsilon(omega);
    %Rh=Rh_fre1(omega);

    if plot_on
        figureParameter
        f1=plot(omega,angle(Rh_fre1(omega,0)),'*r',omega,angle(Rh_fre1(omega,1)),'-og',omega,angle(Rh_fre1(omega,2)),'-.k',omega,-angle(Rx),'-b');
        a1=xlabel('$\omega$');
        xlim([0 1]);
        ylim([0,1.57]);
        set(gca,'Ytick',[0 1.57]);
        a2=ylabel('$\phi_a^*$');
        h1=legend('$B=0$','$B=1$','$B=2$');%,'$\epsilon=0.8$');
        %h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
        set(h1,'location','northeast')
        fig_name='./figure/omega_angleR.eps';
        figurePostTreat
    end

    
    case 2 
        disp("nonlinear response theory:   nonlinear (negative-feedback) adaptive circuit + linear signal")

        Rh_fre1=@(sigma) 1./(1-sqrt(-1)*sigma*tau_h);
        %Rx_fre=@(sigma) sqrt(-1).*R_fre_velo(sigma)./sigma;

        %%%%% this is for negative feedback circuit
        %Q=1;w0=Q;
        %Rx_fre=@(sigma,A) 1./(1+3*A.^2/4-sqrt(-1)*Q*(sigma/w0-w0./sigma)); 

        %% another version include a finite epsilon
        Rx_fre_epsilon=@(sigma,A) 1/gamma_x*1./(1+3*A.^2/4-sqrt(-1)*sigma*tau_a-1./(sqrt(-1)*sigma*tau_y-epsilon));
        %%%%%

        G=@(sigma,A)  Rh_fre1(sigma).*Rx_fre_epsilon(sigma,A);

        %%  plots
        omega=0.00001:0.001:10;
        %omega=0.01.*1.1.^(0:100);

        if plot_on
            figureParameter
            f1=plot(omega,-angle(Rx_fre_epsilon(omega,0)),'-r',omega,-angle(Rx_fre_epsilon(omega,1)),'-og',omega,-angle(Rx_fre_epsilon(omega,2)),'-.k',omega,angle(Rh_fre1(omega)),'*b');
            a1=xlabel('$\omega$');
            xlim([0 1]);
            ylim([0,1.57]);
            set(gca,'Ytick',[0 1.57]);
            %a2=ylabel('$\phi_a^*$');
            h1=legend('$A=0$','$A=1$','$A=2$');%,'$\epsilon=0.8$');
            %h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
            set(h1,'location','northeast')
            fig_name='./figure/omega_angleR.eps';
            figurePostTreat
        end

    case 3 
        disp("nonlinear response theory:   nonlinear (incoherent feedforward) adaptive circuit + linear signal")
        Rh_fre1=@(sigma) 1./(1-sqrt(-1)*sigma*tau_h);
        %Rx_fre=@(sigma) sqrt(-1).*R_fre_velo(sigma)./sigma;

        %%%%% this is for incoherent feedforward circuit
        Rx_fre_epsilon=@(sigma,A) 1/gamma_x*1./(-sqrt(-1)*sigma*tau_a+1+3/4*A.^2)*sqrt(-1).*sigma*tau_y./(sqrt(-1)*sigma*tau_y-1);
        %%%%%

        G=@(sigma,A)  Rh_fre1(sigma).*Rx_fre_epsilon(sigma,A);

        %%  plots
        omega=0.00001:0.001:10;
        %omega=0.01.*1.1.^(0:100);

        if plot_on 
            figureParameter
            f1=plot(omega,-angle(Rx_fre_epsilon(omega,0)),'-r',omega,-angle(Rx_fre_epsilon(omega,1)),'-og',omega,-angle(Rx_fre_epsilon(omega,2)),'-.k',omega,angle(Rh_fre1(omega)),'*b');
            a1=xlabel('$\omega$');
            xlim([0 1]);
            ylim([0,1.57]);
            set(gca,'Ytick',[0 1.57]);
            %a2=ylabel('$\phi_a^*$');
            h1=legend('$A=0$','$A=1$','$A=2$');%,'$\epsilon=0.8$');
            %h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
            set(h1,'location','northeast')
            fig_name='./figure/omega_angleR.eps';
            figurePostTreat
        end
        
end

%%  obtain the oscillation amplitude and frequency from theory


%Pert_array=0:0.1:3; %regime of amplitude perturbation
% the meaning of A and B could be opposite
predicted_B=zeros(length(Pert_array),1);
predicted_A=zeros(length(Pert_array),1);
predicted_w=zeros(length(Pert_array),1);
predicted_N=zeros(length(Pert_array),1);


for j=1:length(Pert_array)
   w_temp=w_ini;
   while -angle(G(w_temp,Pert_array(j)))>0
       w_temp=w_temp+dw;
   end
    
    predicted_w(j)=w_temp;
    predicted_N(j)=1/abs(G(w_temp,Pert_array(j)));
    
    if which_model==1      
    %% linear adaptive circuit
    predicted_B(j)=Pert_array(j);
    predicted_A(j)=predicted_B(j)*alpha2*abs(Rx_fre_epsilon(w_temp)); % A and B are the absolute amplitude of oscillation for a and s, respectively
    else
    %% nonlinear adaptive circuit
    predicted_A(j)=Pert_array(j);
    predicted_B(j)=predicted_A(j)/(alpha2*abs(Rx_fre_epsilon(w_temp,Pert_array(j))));
    end
   
end


% figure, plot(predicted_N,predicted_w,'*r');
% xlim([0 7]);
% ylabel('w');
% 
% figure, plot(predicted_N,predicted_B,'-r');
% xlim([0 7]);
% ylabel('B');
% 
% figure, plot(predicted_B,predicted_w,'-r');
% xlim([0 10]);
% xlabel('B');
% ylabel('w');







