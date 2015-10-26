function [D]=Dcoeff(x,y,N,param_vec)
D=0;
for k=1:N 
   D=D+cos(k*pi*x/1).*cos(k*pi*y/1).*param_vec(k);
end
D=D./N+1+x^2+y^2;   
%D=D+1.7+x^2+y^2; 