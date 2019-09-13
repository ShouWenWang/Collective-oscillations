%% predicted frequency from linear response, varying epsilon
fresh_compute=0;

if fresh_compute
    disp("Fresh computation");
epsilon=0.001*1.2.^(0:40);
max_M=length(epsilon);
tau_star=zeros(1,max_M); % the maximum tau_s that enables collective oscillations
ini_w=0.001; %initial search frequency, parameter used in interpolate the onset frequency

j0=1;
for k=1:max_M
tau_s=cat(2,0.1:0.1:2,2:1:100);
tau_s=sort(tau_s,'descend');
which_model=2;
Pert_array=0; %regime of amplitude perturbation
plot_on=1;
max_N=length(tau_s);
 disp("a new epsilon")

    for j=j0:max_N
    tau_a=1;
    tau_y=1;
    gamma_a=1;
    alpha1=0.5;
    alpha2=0.5;
    para_array=[tau_s(j),tau_a,tau_y,gamma_a,epsilon(k),alpha1,alpha2];
    [predicted_N,predicted_A,predicted_B,predicted_w]=predicting_amp_fre_from_analytical_response_spectrum_adaptation(which_model,para_array,Pert_array,plot_on);
        if predicted_w>ini_w
            tau_star(k)=tau_s(j);
            j0=j;
            disp("found tau_star")
            break;

        end

    end
end

%% save_all_results

tt=today;
aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
output_name="./data/"+aa+"_epsilon_tau_star"+num2str(floor(10*rand()));
save(output_name);

else
    load ./data/190410_epsilon_tau_star_paper

end


index=tau_star>0;
figureParameter
f1=plot(1./epsilon(index),tau_star(index),'*b');
xlim([0 100]);
%ylim([0.6 0.9]);
%ylim([0.2 0.6]);
%set(gca,'YTICK',[0.6 0.7 0.8 0.9]);
%set(gca,'YTICK',[0 0.4 0.8]);
a1=xlabel('$1/\epsilon$');
a2=ylabel('$\tau_s^*$');
fig_name='./figure/tau_epsilon_predict.eps';
figurePostTreat

