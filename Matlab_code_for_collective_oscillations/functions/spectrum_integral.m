function sum_fre=spectrum_integral(spectrum,fre,max_fre)

sum_fre=0;
Max_N=sum(fre<max_fre);
for j=2:Max_N
  sum_fre=sum_fre+(spectrum(j-1)+spectrum(j))/2*(fre(j)-fre(j-1))/pi;
end
