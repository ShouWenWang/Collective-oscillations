
fresh_computation=0;
if fresh_computation
    disp("fresh computation")
%% predicted frequency from linear response, varying tau_s
tau_s=0.01*1.3.^(0:30);


ini_w=0.001; %initial frequency used to infer the intersection 
which_model=2;
Pert_array=0; %regime of amplitude perturbation
plot_on=1;
Max_N=length(tau_s);

predicted_N=zeros(1,Max_N);
predicted_A=zeros(1,Max_N);
predicted_B=zeros(1,Max_N);
predicted_w=zeros(1,Max_N);

for j=1:length(tau_s)
tau_a=1;
tau_y=1;
gamma_a=1;
epsilon=0.1;
alpha1=0.5;
alpha2=0.5;
para_array=[tau_s(j),tau_a,tau_y,gamma_a,epsilon,alpha1,alpha2];
[predicted_N(j),predicted_A(j),predicted_B(j),predicted_w(j)]=predicting_amp_fre_from_analytical_response_spectrum_adaptation(which_model,para_array,Pert_array,plot_on);
end


%% save_all_results

tt=today;
aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/"+aa+"_taus_fre"+num2str(floor(10*rand()));
save(output_name);

else
    load ./data/190410_taus_fre_paper
end


index=predicted_w>ini_w; % there is change!
figureParameter
f1=semilogx(tau_s(index),predicted_w(index),'*b');
xlim([0.01 100]);
ylim([0 1.5]);
%ylim([0.2 0.6]);
%set(gca,'YTICK',[0.6 0.7 0.8 0.9]);
%set(gca,'YTICK',[0 0.4 0.8]);
a1=xlabel('$\tau_s$');
a2=ylabel('$\omega_o$');
fig_name='./figure/tau_w_predict.eps';
figurePostTreat

figureParameter
f1=semilogx(tau_s(index),predicted_N(index),'*b');
xlim([0.01 100]);
ylim([0 20]);
%ylim([0.2 0.6]);
%set(gca,'YTICK',[0.6 0.7 0.8 0.9]);
%set(gca,'YTICK',[0 0.4 0.8]);
a1=xlabel('$\tau_s$');
a2=ylabel('$\bar{N}_o$');
fig_name='./figure/tau_N_predict.eps';
figurePostTreat


