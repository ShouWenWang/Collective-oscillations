function [fre_os_esti,Neff_os_esti]=predicting_amp_fre_from_computed_response_spectrum(Hz,h0,plot_on)

%% second part: estimating the frequency and amplitude based on the response spectrum computed numerically
% the strategy of estimating the frequency and amplitude  is to first
% compute the response spectrum at different perturbation amplitude h0.
% Then, based on Eq. (3) in our main text and using linear interpolation, to infer the precise frequency and
% amplitude at which the collective oscillation. 

fresh_computation=1;

if fresh_computation
    disp("Running fresh computation!!!")
%%%%%%%% deciding the frequency for perturbation
% for the inferrence to be precise, we only need well-resolved spectrum
% around the frequency where oscillations actually happens [this can be
% obtained by collective simulations]
    %Hz=0.03:0.002:0.05; % used in the paper
    %Hz=0.03:0.01:0.05; % for test
 
  
%%%%%% deciding the amplitude of the perturbation
    %h0=0.2:0.2:3; %used in the paper
    %h0=[0.2 0.4]; % for test
    

%%% initialization
M0=length(Hz);
Rx_fre_array=zeros(M0,1);
Rv_fre_array=zeros(M0,1);
Mean_Energy_outflow=zeros(M0,1);
Mean_efficiency=zeros(M0,1);


%wh=1;
%% nonlinear signal
%Rh_fre1=@(sigma,B) 1./(1+3/4*B^2-sqrt(-1)*sigma); % the signal spectrum, assuming a qubic nonlinear term 
%% linear signal
Rh_fre1=@(sigma,B) 1./(1-sqrt(-1)*sigma); % linear version 


which_model=2;
for l=1:length(h0)
 for j=1:length(Hz)%M0
   [Rx_fre_array(j),Rv_fre_array(j),time,h,x_bar,Energy_outflow,efficiency,Internal_power]=Computing_response_from_direct_simulation(Hz(j),h0(l),which_model,plot_on);
   Mean_Energy_outflow(j)=mean(Energy_outflow);
   Mean_efficiency(j)=mean(efficiency);
 end
 

    Neff=1./abs(Rx_fre_array.*Rh_fre1(2*pi*Hz,h0(l))');    
    data_energyout(l,:)={-angle(Rx_fre_array),Rx_fre_array,Rv_fre_array,Mean_Energy_outflow, Mean_efficiency,Neff}; 
end



 
fre0=2*pi*Hz;
fre_os_esti=zeros(1,length(h0));
Neff_os_esti=zeros(1,length(h0));

for l=1:length(h0)
    angle1=data_energyout{l,1}-angle(Rh_fre1(2*pi*Hz,h0(l)))';
    Neff=data_energyout{l,6};

    j=1;
    while angle1(j)>0
      j=j+1;
    end
    
    k0=(angle1(j-1)-angle1(j))/(fre0(j-1)-fre0(j));
    fre_os_esti(l)=fre0(j)-angle1(j)/k0;
    
    k1=(Neff(j-1)-Neff(j))/(fre0(j-1)-fre0(j));
    Neff_os_esti(l)=Neff(j)+k1*(fre_os_esti(l)-fre0(j));
    
end


%% save_all_results

% tt=today;
% aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/predicted_amp_fre_"+num2str(floor(10*rand()));
save(output_name);




else
    
    load  ./data/190406_predicted_amp_fre_1
    
end


if plot_on


if length(h0)>=3
    
    % index of the perturbation amplitude to look at
    index=[1 2 3]; % default choice

else if length(h0)==2
        
        index=[1 1 2];
    
    else
        index=[1 1 1];
    end
end
    

angle_x1=data_energyout{index(1),1};
angle_x2=data_energyout{index(2),1};
angle_x3=data_energyout{index(3),1};

h0_1=h0(index(1));
h0_2=h0(index(2));
h0_3=h0(index(3));

Omega=0:0.001:10;
angle_h1=angle(Rh_fre1(Omega,h0_1));
angle_h2=angle(Rh_fre1(Omega,h0_2));
angle_h3=angle(Rh_fre1(Omega,h0_3));

% G=1/|RxRh|
G1=1./abs(data_energyout{index(1),2}.*Rh_fre1(2*pi*Hz,h0_1)');
G2=1./abs(data_energyout{index(2),2}.*Rh_fre1(2*pi*Hz,h0_2)');
G3=1./abs(data_energyout{index(3),2}.*Rh_fre1(2*pi*Hz,h0_3)');



figureParameter
f1=plot(Hz*2*pi,angle_x1,'-*r',Hz*2*pi,angle_x2,'-.k',Hz*2*pi,angle_x3,'-og');%,Hz*2*pi,angle4,'-+b',Hz*2*pi,angle5,'-.cyan',Omega,angle(Rh_fre1(Omega)),'-k');
xlim([0 0.4]);
ylim([0 1.6]);
set(gca,'YTICK',[0 1.57]);
a1=xlabel('$\omega$');
a2=ylabel('$\phi_a^*$');
h1=legend(strcat("$B=",num2str(h0_1),"$"), strcat("$B=",num2str(h0_2),"$"),strcat("$B=",num2str(h0_3),"$"));
fig_name='./figure/phase_shift.eps';
legend boxoff
figurePostTreat



figureParameter
f1=plot(Hz*2*pi,angle_x1,'-*r',Hz*2*pi,angle_x2,'-.k',Hz*2*pi,angle_x3,'-og',...
    Omega,angle_h1,'--r',Omega,angle_h2,'--k',Omega,angle_h3,'--g');
xlim([0 0.4]);
ylim([0 1.6]);
set(gca,'YTICK',[0 1.57]);
a1=xlabel('$\omega$');
a2=ylabel('$\phi_a^*$');
h1=legend(strcat("$\phi_a: B=",num2str(h0_1),"$"), strcat("$\phi_a:B=",num2str(h0_2),"$"),strcat("$\phi_a:B=",num2str(h0_3),"$"),strcat("$\phi_s: B=",num2str(h0_1),"$"), strcat("$\phi_s: B=",num2str(h0_2),"$"),strcat("$\phi_s: B=",num2str(h0_3),"$"));
fig_name='./figure/angle_xh_intersection.eps';
legend boxoff
figurePostTreat



figureParameter
f1=plot(Hz*2*pi,G1,'-*r',Hz*2*pi,G2,'-.k',Hz*2*pi,G3,'-og');
xlim([0 0.4]);
ylim([0 1.6]);
set(gca,'YTICK',[0 1.57]);
a1=xlabel('$\omega$');
a2=ylabel('$1/|\tilde{R}_a\tilde{R}_h|$');
h1=legend(strcat("$B=",num2str(h0_1),"$"), strcat("$B=",num2str(h0_2),"$"),strcat("$B=",num2str(h0_3),"$"));
fig_name='./figure/effective_N.eps';
legend boxoff
figurePostTreat



% figureParameter
% f1=semilogx(h0,Mean_efficiency,'-or');
% a1=xlabel('$h_0$');
% xlim([0.1 100]);
% set(gca,'XTICK',[0.1 1 10 100]);
% a2=ylabel('$\eta$');
% %a2=ylabel('$\eta$');
% %a2=ylabel('$-\dot{W}_{ext}$');
% %h1=legend('Exact','$h_0=0.5$','$h_0=1$','$h_0=5$');
% %set(h1,'location','northwest');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% figureParameter
% f1=plot(h0,Mean_efficiency,'-or');
% %a1=xlabel('$h_0$');
% xlim([0 10]);
% %a2=ylabel('$\dot{W}_{int}$');
% %a2=ylabel('$\eta$');
% %a2=ylabel('$-\dot{W}_{ext}$');
% %h1=legend('Exact','$h_0=0.5$','$h_0=1$','$h_0=5$');
% %set(h1,'location','northwest');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% figureParameter
% f1=semilogx(h0,Rv_fre_array','-or');
% a1=xlabel('$h_0$');
% xlim([0.01 100]);
% a2=ylabel('$\tilde{R}_v''$');
% %h1=legend('Exact','$h_0=0.5$','$h_0=1$','$h_0=5$');
% %set(h1,'location','northwest');
% fig_name='./figure/efficiency.eps';
% figurePostTreat

% 
% 
% %% test whether the numerically simulated response spectrum is correct 


% k=1;
% Rx_fre_array_1=data_energyout{k,1};
% Rv_fre_array_1=data_energyout{k,2};
% Mean_Energy_outflow_1=data_energyout{k,3};
% Mean_efficiency_1=data_energyout{k,4};
% 
% k=2;
% Rx_fre_array_2=data_energyout{k,1};
% Rv_fre_array_2=data_energyout{k,2};
% Mean_Energy_outflow_2=data_energyout{k,3};
% Mean_efficiency_2=data_energyout{k,4};
% 
% k=3;
% Rx_fre_array_3=data_energyout{k,1};
% Rv_fre_array_3=data_energyout{k,2};
% Mean_Energy_outflow_3=data_energyout{k,3};
% Mean_efficiency_3=data_energyout{k,4};



% Omega1=2*pi*Hz;
% 
% figureParameter
% f1=semilogx(Omega,real(data_cell{1,2})./Omega,'-k',Omega1,Rx_fre_array_1,'*r',Omega1,Rx_fre_array_3,'+b',Omega1,Rx_fre_array_2,'og');
% a1=xlabel('$\omega$');
% xlim([0.1 10]);
% a2=ylabel('$\tilde{R}_x''''$');
% h1=legend('Exact','$h_0=0.5$','$h_0=1$','$h_0=5$');
% set(h1,'location','northwest');
% figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega,real(data_cell{1,2}),'-k',Omega1,Rv_fre_array_1,'*r',Omega1,Rv_fre_array_3,'+b',Omega1,Rv_fre_array_2,'og');
% a1=xlabel('$\omega$');
% xlim([0.01 100]);
% a2=ylabel('$\tilde{R}_v''$');
% h1=legend('Exact','$h_0=0.5$','$h_0=1$','$h_0=5$');
% set(h1,'location','northwest');
% figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega1,Mean_efficiency_1/(0.5^2),'-or',Omega1,Mean_efficiency_3/(1^2),'-+b',Omega1,Mean_efficiency_2/(5^2),'-*g');
% xlim([0.1 10]);
% a1=xlabel('$\omega$');
% a2=ylabel('$\eta/h_0^2$');
% h1=legend('$h_0=0.5$','$h_0=1$','$h_0=5$');
% set(h1,'location','northeast');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega1,Mean_Energy_outflow_1./Mean_efficiency_1,'-or',Omega1,Mean_Energy_outflow_3./Mean_efficiency_3,'-+b',Omega1,Mean_Energy_outflow_2./Mean_efficiency_2,'-*g');
% xlim([0.1 10]);
% ylim([5  15]);
% a1=xlabel('$\omega$');
% a2=ylabel('$\dot{W}_{int}$');
% h1=legend('$h_0=0.5$','$h_0=1$','$h_0=5$');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega1,Mean_Energy_outflow_1/0.5^2,'-or',Omega1,Mean_Energy_outflow_3,'-+b',Omega1,Mean_Energy_outflow_2/5^2,'-*g');
% xlim([0.1 10]);
% a1=xlabel('$\omega$');
% a2=ylabel('$-\dot{W}_{ext}/h_0^2$');
% h1=legend('$h_0=0.5$','$h_0=1$','$h_0=5$');
% set(h1,'location','northeast');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% figureParameter
% f1=plot(Mean_Energy_outflow_1/0.5^2,Mean_efficiency_1/0.5^2,'-or',Mean_Energy_outflow_3,Mean_efficiency_3,'-+b',Mean_Energy_outflow_2/5^2,Mean_efficiency_2/5^2,'-*g');
% xlim([-0.5 0.5]);
% a1=xlabel('$-\dot{W}_{ext}/h_0^2$');
% a2=ylabel('$\eta/h_0^2$');
% h1=legend('$h_0=0.5$','$h_0=1$','$h_0=5$');
% set(h1,'location','southeast');
% fig_name='./figure/efficiency.eps';
% figurePostTreat
% 
% 
% % figureParameter
% % f1=plot(-10:0.1:20,-10:0.1:20,'-b',h0^2*Rv_fre_array/2,-Mean_Energy_outflow,'or');
% % a1=xlabel('$h_0^2\omega\tilde{R}_x''''/2$');
% % xlim([-5 15]);
% % ylim([-5 15]);
% % a2=ylabel('$\dot{W}_{ext}$');
% % %h1=legend('$\omega\tilde{C}_x$','$2\tilde{R}_x''''$');
% % fig_name='./figure/response_compare_power1.eps';
% % figurePostTreat
% 
% % figureParameter
% % f1=semilogx(Omega1,h0^2*Rv_fre_array/2,'-b',Omega1,-Mean_Energy_outflow,'or');
% % a1=xlabel('$\omega$');
% % xlim([0.01  100]);
% % h1=legend('$h_0^2\omega\tilde{R}_x''''/2$','$\dot{W}_{ext}$');
% % set(h1,'location','northwest');
% % fig_name='./figure/response_compare_power2.eps';
% % figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega, real(data_cell{1,2})./Omega,'-r',Omega1,Rx_fre_array,'-*b');
% a1=xlabel('$\omega$');
% xlim([0.01 100]);
% a2=ylabel('$\tilde{R}_x''''$');
% h1=legend('$h_0=0.5$','$h_0=5$');
% set(h1,'location','northwest');
% fig_name='./figure/response_compare_exact_simu.eps';
% figurePostTreat
% 
% figureParameter
% f1=semilogx(Omega, real(data_cell{1,2}),'-r',Omega1,Rv_fre_array,'-*b');
% a1=xlabel('$\omega$');
% xlim([0.01 100]);
% a2=ylabel('$\tilde{R}_v''$');
% h1=legend('$h_0=0.5$','$h_0=5$');
% set(h1,'location','northwest');
% fig_name='./figure/response_velo_compare_exact_simu.eps';
% figurePostTreat

% figure,semilogx();
% xlim([0.1 100]);

end