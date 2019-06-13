function [Omega,Rx_fre,Rv_fre]=nonlinear_response(h0,which_model)

% h0 is the oscillation amplitude
% the algorithm calculates the response for all different frequencies
% which_model=1;
%% common simulation parameters  

dt=0.01;
T0=500000;

L=floor(T0/dt);

time=(1:L)*dt;
%Hz=0.026;
Hz=0.01.*floor(1.3.^(1:20));%

%Hz=0.1.*floor(1.3.^(1:20));%

h=0;
for i=1:length(Hz)
h=h+h0*cos(2*pi*Hz(i)*time);
end

switch which_model
    case 1  %%  adaptive model
        disp("Running adaptive model")
        wx=1;
        wy=1;
        epsilon=0.1;
        %wyy=0.1;
        gamma_x=1;
        gamma_y=1;
        T=1;


        x=zeros(L,1);
        y=zeros(L,1);

        y(1)=2; %initial condition
        x(1)=1;

        noise_x=randn(L,1)*sqrt(dt*2*T/gamma_x);
        noise_y=randn(L,1)*sqrt(dt*2*T/gamma_y);



        for j=2:L

            x(j)=x(j-1)-wx*(x(j-1)-y(j-1))*dt+h(j-1)*dt+noise_x(j-1)-x(j-1)^3*dt;
            y(j)=y(j-1)-wy*x(j-1)*dt-wy*epsilon*y(j-1)*dt+noise_y(j-1);
        end

    case 2  %%  dicty/cAMP
        disp("Running FHN model")
        wy=0.2;
        epsilon=0.1;
        x0=1.5;
        %self-degradation rate

        gamma_x=1;
        gamma_y=1;
        T=0.1;

        x=zeros(L,1);
        y=zeros(L,1);



        noise_x=randn(L,1)*sqrt(dt*2*T/gamma_x);
        noise_y=randn(L,1)*sqrt(dt*2*T/gamma_y);


        for j=2:L

            x(j)=x(j-1)+(x(j-1)-(x(j-1)).^3/3)*dt-y(j-1)*dt+noise_x(j-1)+h(j-1)*dt;
            y(j)=y(j-1)+wy.*(x(j-1)-epsilon.*y(j-1)+x0)*dt+noise_y(j-1);
        end

    case 3  %% ion polarity
        disp("Running ion polarity model")
        KB=5*exp(-2*h);
        KF=5*exp(2*h);

        k41=250;
        k14=50;
        k32=10;
        k23=2;
        SL=1;
        SR=10;

        Pa=zeros(L,1);
        Pb=zeros(L,1);
        x=zeros(L,1);
        Pa(1)=0.5;
        Pb(1)=0.5;

        for i=1:L-1
            P1=Pa(i)*KB(i)/(KB(i)+1);
            P2=Pa(i)/(KB(i)+1);
            P3=Pb(i)*KF(i)/(KF(i)+1);
            P4=Pb(i)/(KF(i)+1);
            x(i)=P2+P3-P1-P4;

            Pa(i+1)=Pa(i)+dt*(P4*k41*SL+P3*k32-P1*k14-P2*k23*SR);
            Pb(i+1)=Pb(i)+dt*(P1*k14+P2*k23*SR-P4*k41*SL-P3*k32);

        end

end

%% figure

% figureParameter
% f1=plot(time,h,'-r',time,x,'--b',time,y,'-.k');
% a1=xlabel('time');
% xlim([20 100]);
% %ylim([-10 15]);
% h1=legend('$h$','$x$','$y$');
% set(h1,'location','southwest');
% fig_name='./figure/Adaptation-oscillation-noise.eps';
% figurePostTreat;

%%  spectrum 

X_fre=(fft(x,L)/L);
h_fre=fft(h,L)/L;
fre=2*pi*(0:L/2)'/(L*dt); 

% figure, semilogx(fre,abs(X_fre(1:L/2+1)),'-r',fre,abs(h_fre(1:L/2+1)),'-.b')
% xlabel('frequency');
% ylim([0.001 1]);
% ylabel('|$\tilde{x}|$','interpreter','latex');
% legend('x','h');

% figureParameter
% f1=semilogx(fre,real(X_fre(1:L/2+1)),'-r',fre,abs(h_fre(1:L/2+1)),'-.b');
% a1=xlabel('$\omega$');
% xlim([0 1]);
% %ylim([0.001 1]);
% a2=ylabel('$\tilde{x}''$');
% fig_name='./figure/Adaptation-nonlinear-x-real.eps';
% figurePostTreat;

% %legend('x','h');
% 
% figureParameter
% f1=semilogx(fre,imag(X_fre(1:L/2+1)),'-r',fre,abs(h_fre(1:L/2+1)),'-.b');
% xlim([0 1]);
% %ylim([0.001 1]);
% a1=xlabel('$\omega$');
% ylim([0.001 1]);
% a2=ylabel('$\tilde{x}''''$');
% fig_name='./figure/Adaptation-nonlinear-x-imag.eps';
% figurePostTreat;


Omega=2*pi*Hz;  %perturbed frequency


index=round(L*dt*Hz+1); % index of given perturbed frequency
Y=X_fre(index);     
Z=h_fre(index);
Rx_fre=Y'./Z;
Rv_fre=-sqrt(-1)*Omega.*Rx_fre;


%% calibration of the aliasing effect % this calibaration is needed only when dt is the simulated dt.
%%  therefore, this step is not needed in experiment
alias_effect=@(x) (2*pi*x).*sin(2*pi*x)./(2-2*cos(2*pi*x)); % x=fre/fre_sampling, range: 0: pi
%dt_simu=0.05*0.01; % this is the update time in simulation
fre_sp=2*pi/dt; %sampling omega
Rx_fre=Rx_fre./alias_effect(2*pi*Hz/fre_sp);
Rv_fre=Rv_fre./alias_effect(2*pi*Hz/fre_sp);


figureParameter
semilogx(Omega,real(Rx_fre),'+r',Omega,imag(Rx_fre),'*b')
a1=xlabel('Frequency');
a2=ylabel('Spectrum');
h1=legend('$\tilde{R}_x$:real','$\tilde{R}_x$:imag');
fig_name='./figure/nonlinear_response.eps';
figurePostTreat
