%% purpose, 19-9-4
% compute the oscillation amplitude an frequency via simulating the population dynamics as the cell density is changed


fresh_computation=1;
if fresh_computation
    %dt=0.01;
    T0=20; %default 20, maximum number of period that you want to use to compute frequency and amplitude
    Time_ini=3; % the transient time, after which we start to compute the amplitude and frequency of the trajectory
    Time_window=0.1; % the sliding window to pick peaks or troughs, this should be around the period of oscillation  
    step_size=0.0001;
    plot_on=1;
    density=cat(2,[0.001:0.02:1,1]); 
    %density=0.001:0.02:0.4; 

    with_explicit_diffusion=1;
    D=0.4;  % 0.4 is the constant at which DQS could end at high cell density
    dt=0.00001; % for D=100, you need to use dt=0.000002 to avoid blowing-up solution

    fre_final_a=zeros(1,length(density));
    amp_final_a=zeros(1,length(density));

    fre_final_s=zeros(1,length(density));
    amp_final_s=zeros(1,length(density));

    ave_level_s=zeros(1,length(density));
    ave_level_s_ex=zeros(1,length(density));
    ave_level_a=zeros(1,length(density));


    plot_fig=1;
    L_mid=floor(2/dt);

    parfor j=1:length(density)

    if with_explicit_diffusion
        [s_auto,s_ex_auto,a_auto,time]=minimum_glycolysis(dt,T0,density(j),D,plot_fig);
    else
        [s_auto,a_auto,time]=minimum_glycolysis_no_Diffusion(dt,T0,density(j),plot_fig);
    end

    ave_level_a(j)=mean(a_auto(L_mid:end));
    ave_level_s(j)=mean(s_auto(L_mid:end));
    ave_level_s_ex(j)=mean(s_ex_auto(L_mid:end));

    [rela_unc_mean_fre,rela_unc_mean_amp,fre_final_s(j),amp_final_s(j)]=oscillation_amplitude_fre_estimation(s_auto,dt,T0,Time_window,plot_fig,Time_ini);
    %[rela_unc_mean_fre,rela_unc_mean_amp,fre_final_a(j),amp_final_a(j)]=oscillation_amplitude_fre_estimation(a_auto,dt,T0,Time_window,plot_fig);


    %[fre_final_s(j),amp_final_s(j)]=oscillation_amplitude_fre_estimation_only_count_peaks(s_auto,dt,T0,Time_window,plot_fig);
    %[fre_final_a(j),amp_final_a(j)]=oscillation_amplitude_fre_estimation_only_count_peaks(a_auto,dt,T0,Time_window,plot_fig);
    disp(num2str(density(j)));
    end
    

    file_name=strcat('./data/glycolysis_amp_fre_density_',num2str(floor(rand(1)*10)));
    save(file_name,'perturb_fre','fre_final_s','amp_final_s');
else
    load ./data/glycolysis_amp_fre_density_1
end

pred_density=0.73;
pred_amp=0;
pred_fre=125;

figureParameter
f1=plot(density,amp_final_s,'-b',pred_density,pred_amp,'or');
ylim([0 0.1])
xlim([0 1])
xlabel("Cell density")
ylabel('Oscillation amplitude');
h=legend('Simulation','Prediction');
set(h,'location','northeast')
fig_name="./figure/density_amp_glycolysis_s.eps";
figurePostTreat

figureParameter
f1=plot(density,fre_final_s,'-b',pred_density,pred_fre,'or');
xlim([0 1])
ylim([100 130])
xlabel("Cell density")
ylabel('Oscillation frequency');
h=legend('Simulation','Prediction');
set(h,'location','southeast')
fig_name="./figure/density_fre_glycolysis_s.eps";
figurePostTreat
