%%  solve for sigma (the instability) with linear response functions, need to provide R_fre_velo, which is generated through exact numerical spectrum
wh=1;tau=0;
Rh_fre=@(sigma) 1./(wh-sqrt(-1)*sigma).*exp(sqrt(-1)*tau*sigma);
Rx_fre=@(sigma) sqrt(-1).*R_fre_velo(sigma)./sigma;
G=@(sigma)  Rh_fre(sigma).*Rx_fre(sigma);

omega_0=0.005*(120:180);
%temp_wr=wr_0; dw=-0.001;
%wr_0=10^(-5);
%temp_wr=temp_wr-dw;
correct_sigma=zeros(1,20);
correct_N=zeros(1,20);
k=1; 
for a=1:length(omega_0)
    temp_wr=omega_0(a);

   
   u=0.001*(-5000:5000);
   sigma=temp_wr+sqrt(-1)*u;
   phi0=-angle(G(sigma));
   for m0=1:length(sigma)-1
       if   phi0(m0)*phi0(m0+1)<0
           correct_sigma(k)=(sigma(m0)+sigma(m0+1))/2;
           correct_N(k)=1/(G(correct_sigma(k)));
           k=k+1;
       end
   end
   
end

index=imag(correct_sigma)>-1;
figureParameter
%f1=plot(real(result),imag(result),'-r',Omega(1:end-1),u1,'-.b');
f1=plot(real(correct_sigma(index)),imag(correct_sigma(index)),'*r');
a1=xlabel('$\omega$');
a2=ylabel('$u$');
fig_name='./figure/omega_u.eps';
figurePostTreat


figureParameter
f1=plot(correct_N(index),imag(correct_sigma(index)),'or');
xlim([-10 10]);
%ylim([-1 1]);
a1=xlabel('$\bar{N}$');
a2=ylabel('$u$');
fig_name='./figure/N_u.eps';
figurePostTreat

figureParameter
f1=plot(correct_N(index),real(correct_sigma(index)),'or');
xlim([-10 10]);
%ylim([-1 1]);
a1=xlabel('$\bar{N}$');
a2=ylabel('$\omega$');
fig_name='./figure/N_w.eps';
figurePostTreat

index=imag(correct_sigma)>0 & imag(correct_sigma)<0.2;



scatter(real(correct_sigma),imag(correct_sigma),'*r');
figure,
xlim([-1 3]);

figure,scatter(correct_N,G(correct_sigma).*correct_N,'*r');