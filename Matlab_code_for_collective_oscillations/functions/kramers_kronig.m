function [Rx_real_cal,Rx_imag_cal,Rx_imag_cal_2]=kramers_kronig(w)
a=1; 
Rx_real=@(x) x*(x<1)+1*(x>=1);
Rx_imag=@(x) x/(a^2+x.^2);
f=@(x) Rx_imag(x).*x./(x.^2-w.^2);
g=@(x) -w*Rx_real(x)./(x.^2-w.^2);

dx=0.01;

L=1000;
Rx_real_cal=0;
Rx_real_theo=a/(a^2+w^2);
Rx_imag_theo=w/(a^2+w^2);
% for i=0:round(w/dx)-1
%     tot=tot+f((i+1/2)*dx)*dx*2/pi;
% end

for i=0:round(L/dx)
    Rx_real_cal=Rx_real_cal+f((i+1/2)*dx)*dx*2/pi;
end

Rx_imag_cal=0;
for i=0:round(L/dx)
    Rx_imag_cal=Rx_imag_cal+g((i+1/2)*dx)*dx*2/pi;
end

Rx_imag_cal_2=0;
for i=round(2*w/dx):round(L/dx)
    Rx_imag_cal_2=Rx_imag_cal_2+g((i+1/2)*dx)*dx*2/pi;
end



%% calculate the curvature
% tot_0=0;
% f0=@(x) x./(a^2+x.^2)./x;
% for i=0:round(L/dx)
%     tot_0=tot_0+f0((i+1/2)*dx)*dx*2/pi;
% end
% 
% tot_1=0;
% f1=@(x) x./(a^2+x.^2)./x;
% for i=0:round(L/dx)
%     tot_1=tot_1+f1((i+1/2)*dx)*dx*2/pi;
% end
% 
% tot_2=0;
% % f2=@(x) x./(a^2+x.^2).*(x./(x.^2-dx.^2)-1./x);
% f2=@(x) x./(a^2+x.^2).*(x./(x.^2-dx.^2)-1./x);
% for i=0:round(L/dx)
%     tot_2=tot_2+f2((i+1/2)*dx)*dx*2/pi;
% end
% 
% 
% tot2=0;
% g=@(x) x^3.*(x./(x.^2-4*dx^2)-2.*x./(x.^2-dx^2)+1./x);
% for i=0:round(1000/dx);
%     tot2=tot2+g((i+1/2)*dx)*dx*2/pi;
% end
% 
% tot2/(dx^2)
% 
% real_R=@(w) a/(a^2+w^2);
% (real_R(2*dx)+real_R(0)-2*real_R(dx))/(dx^2)



% 
% 
% x=1.5:dx:10;
% figure,
% plot(x,g(x),'-r');


