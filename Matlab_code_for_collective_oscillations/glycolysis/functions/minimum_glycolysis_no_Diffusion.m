function [s,a,time]=minimum_glycolysis_no_Diffusion(dt,T0,phi,plot_fig)

%% updated, 19-9-2
% the model is rigid, and a very small dt is needed to ensure accuracy.
% 

% a, PYR
% p, ATP
% y, TRIO etc
% z, BPG etc
% s, Ace
% phi, cell density
%%


%% parameters
h=3; alpha1=1; alpha2=1; epsilon=0.01; Kin=0.5; Kex=0.3; tau=0.01; c0=0.02; tau_s=0.001;



% dt=0.00001;
% T0=10;
L=floor(T0/dt);
time=dt*(0:L-1)';

v1=@(p) p/(1+p^(2*h));
v2=@(p,z) z/(1+p^2);
v3=@(p) p^2/(1+p^2);

y=zeros(L,1);
z=zeros(L,1);
p=zeros(L,1);
a=zeros(L,1);
s=zeros(L,1);


s0=phi/(phi*Kin+Kex);
p(1,1)=1; % initial condition
z(1,1)=1;
y(1,1)=1/(alpha2*s0+c0); % initial condition
a(1,1)=1/alpha1;
s(1,1)=s0;


for i=1:L-1
    y(i+1)=y(i)+dt/tau*(2*v1(p(i))-(alpha2*s(i)+c0)*y(i)-epsilon*y(i));
    z(i+1)=z(i)+dt/tau*((alpha2*s(i)+c0)*y(i)-2*v2(p(i),z(i)));
    p(i+1)=p(i)+dt/tau*(-2*v1(p(i))+4*v2(p(i),z(i))-2*v3(p(i)));
    a(i+1)=a(i)+dt/tau*(2*v2(p(i),z(i))-alpha1*a(i));
    s(i+1)=s(i)+dt/tau_s*(phi/(1+phi)*alpha1*a(i)-(phi*Kin+Kex)/(1+phi)*s(i));



end


if plot_fig    
    
    
    %% DIY colors for the plots
    red = [255 0 0]/255;
    dark_green = [0 180 0]/255;
    light_green=[0 255 0]/255;
    blue=[0 0 255]/255;
    grey=[150 150 150]/255;
    orange=[255 178 102]/255;
    dark_blue='#0000CC';
    light_red='#FF6666';
    light_blue='#6666FF';


    figureParameter
    %f1=plot(time,s,'-r',time,a-1,'-b',time,(y-1/s0)*s0,'-g',time,z-1,'-k',time,p-1,'-cyan');
    %plot(time,s_in,'-r',time,s_ex,'-.m',time,p,'-b',time,y,'-g',time,z,'-k',time,a,'-cyan');
    f1=plot(time,s,'-',time,p,'--',time,y,'-.',time,z,'-',time,a,'-');
    xlabel('Time');
    %xlim([9 10]);
    %ylim([-0.2, 0.2]);
    set(f1,{'color'},{blue;grey;dark_green;orange;red})


    %h1=legend('$s$','$p$','$y$','$z$','$a$');
%     %h1=legend('$\delta s$','$\delta a$','$\delta y$','$\delta z$','$\delta o$');%'$x_1$','$x_2$');%,'$x_2$','$y_2$');
%     %h1=legend('S4int','S4ext','N2');
    %set(h1,'location','best');
    fig_name='./figure/ACA_NADH.eps';
    figurePostTreat;


    % figureParameter
    % f1=plot(time,s-5,'-r');
    % a1=xlabel('time');
    % xlim([0 T0]);
    % ylim([-0.01, 0.01]);
    % %h1=legend('$s$','$a$','$y$','$z$','$o$');%'$x_1$','$x_2$');%,'$x_2$','$y_2$');
    % %h1=legend('S4int','S4ext','N2');
    % %set(h1,'location','northwest');
    % fig_name='./figure/ACA_NADH.eps';
    % figurePostTreat;

    % T_min=0.2;
    % [rela_unc_mean_fre,rela_unc_mean_amp,fre_final,amp_final]=oscillation_amplitude_fre_estimation(s,dt,T0,T_min);
    
end
