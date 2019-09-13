%% purpose, 19-9-4
% compute the response spectrum of the cell to external signal, 
% and thus computationally determine the onset/end condition via self-consistency

startup

clear all; close all;

fresh_computation=1;
if fresh_computation

    %% numerical results
    %dt=0.01;
    max_period=100; %default 20, maximum number of period that you want to use to compute frequency and amplitude
    step_size=0.0001;
    glycolytic_model=1; % 1, for the glycolytic circuit; 0 for the adaptive system, used as a control
    plot_fig=1;
    density=0.73; 
    D=0.4;
    
    plot_fig_0=1;
    [s_in_auto,s_ex_auto,time]=minimum_glycolysis(0.00001,10,density,D,plot_fig_0);
    s_ex_0=mean(s_ex_auto); %0.0737; % the static level of the external signal, determined from the autonomous simulation


    if glycolytic_model
        perturb_fre=cat(2,[1:5:120,120:1:130,132:5:200]);
%        perturb_fre=122:0.5:127;
    
    else
        perturb_fre=0.01*1.2.^(0:30);
    end

    tot_N=length(perturb_fre);

    phase_shift=zeros(tot_N,1);
    Response_amp=zeros(tot_N,1);
    dt=0; % just to initialize this variable

    parfor j=1:length(perturb_fre)

    T0=2*pi*max_period/perturb_fre(j);
    dt=0.00001;%step_size*2*pi/perturb_fre(j);

    if glycolytic_model
        %[x,h,time]=synthetic_circuit_periodic_perturbation(dt,T0,perturb_fre(j),density,He_static,plot_fig);
        [x,h,time]=minimum_glycolysis_periodic_perturbation(dt,T0,perturb_fre(j),D,s_ex_0,plot_fig);

    else
        [x,h,time]=adaptation_perturbation(dt,T0,perturb_fre(j));
    end

    [phase_shift(j),x_amp,h_amp]=detect_phase_shift(perturb_fre(j),h,dt,x);
    Response_amp(j)=x_amp/h_amp;

    end

    xx=num2str(dt);
    yy=num2str(density);
    file_name=strcat('./data/glycolysis_response_amp_phi_MaxPeriod',num2str(max_period),'_dt',xx(3:end),'_d',yy(3:end));
    save(file_name,'perturb_fre','Response_amp','phase_shift','density');
else
    load ./Data/glycolysis_response_amp_phi_MaxPeriod100_dt_d73
    D=0.4;
end



figureParameter
f1=plot(perturb_fre,phase_shift,'*b');%,perturb_fre,os_fre_x,'-*b',perturb_fre,os_fre_h2,'-*g');
a1=xlabel('$\omega$');
ylim([-1.6, 1.6]);
% xlim([0 0.3]);
a2=ylabel('Phase shift of $s_{in}$');
set(gca,'Ytick',[-1.6 0 1.6]);
%set(gca,'Xtick',[0.01 0.1 1]);
%h1=legend('New method:$h\;\;$','New method:$x\;\;$','Old method:$h\;$');
%set(h1,'location','southeast')
%ylim([0 0.3]);
fig_name='./figure/phase_shift.eps';
figurePostTreat

figureParameter
f1=plot(perturb_fre,Response_amp,'*b');%,perturb_fre,os_fre_x,'-*b',perturb_fre,os_fre_h2,'-*g');
a1=xlabel('$\omega$');
ylim([0, 2]);
% xlim([0 0.3]);
a2=ylabel('Response amplitude of $s_{in}$');
%set(gca,'Xtick',[0.01 0.1 1]);
%h1=legend('New method:$h\;\;$','New method:$x\;\;$','Old method:$h\;$');
%set(h1,'location','southeast')
%ylim([0 0.3]);
fig_name='./figure/response_amp.eps';
figurePostTreat

%ylim([0 0.3]);
%xlim([0 3]);

%% for the glycolysis
tau_s=0.0001;
Kex=0.3;
Rh_fre1=@(omega) density*D./(Kex-sqrt(-1)*omega*tau_s+density*D);%.*exp(sqrt(-1)*tau*sigma);

%omega=0.001:0.02:1;
omega=0.01.*1.2.^(0:100);
%Rh=Rh_fre1(omega);

figureParameter
f1=plot(perturb_fre,angle(Rh_fre1(perturb_fre)),'-b',perturb_fre,phase_shift,'-r');
%f1=semilogx(omega,angle(Rh_fre1(omega)),'-.g',perturb_fre,phase_shift,'-r');
xlabel('Frequency');
%xlim([0 0.45]);
ylim([-1.58,1.58]);
%set(gca,'Ytick',[-1.57 0 1.57]);
yticks([-pi/2  0 pi/2])
yticklabels({'-\pi/2','0','\pi/2'})
ylabel('Phase shift');
h1=legend('-$\phi_{s_{ex}}$','$\phi_{s_{in}}$');%,'$\epsilon=0.8$');
%h1=legend('$Q=0.5$','$Q=1$','$Q=2$','$Q=3$');
set(h1,'location','northwest')
fig_name='./figure/omega_angleR.eps';
figurePostTreat

compound_response_amp=Response_amp'.*abs(Rh_fre1(perturb_fre));
target=ones(1,length(compound_response_amp));
figureParameter
%plot(perturb_fre,compound_response_amp,'-g',perturb_fre,target,'--r');
f1=plot(perturb_fre,compound_response_amp,'-r');
xlabel('Frequency');
%set(gca,'Ytick',[-1.57 0 1.57]);
a2=ylabel('$|\tilde{R}_{s_{in}}\tilde{R}_{s_{ex}}|$');
%h=legend('$|\tilde{R}_a\tilde{R}_s|$','');
fig_name='./figure/omega_compound_amp.eps';
figurePostTreat

%% determine the intersection


% higher end
line1_x1=124.5;
line1_y1=1.0464;
line1_x2=125;
line1_y2=0.98237;

line2_x1=124.5;
line2_y1=1;
line2_x2=125;
line2_y2=1;


k1=(line1_y2-line1_y1)/(line1_x2-line1_x1);



k2=(line2_y2-line2_y1)/(line2_x2-line2_x1);

%%% the correct frequency
solution=((line2_y1-line1_y1)-(k2*line2_x1-k1*line1_x1))/(k1-k2);


