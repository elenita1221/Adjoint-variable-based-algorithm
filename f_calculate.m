function [f]=f_calculate
syms x y
u=0;N=10;param_vec(1)=3;param_vec(2)=5;param_vec(3)=7;param_vec(4)=11;param_vec(5)=13;
param_vec(6)=17;param_vec(7)=19;param_vec(8)=23;param_vec(9)=29;param_vec(10)=31;
for k=1:N 
   u=u+sin(k*pi*x/1).*sin(k*pi*y/1).*param_vec(k);
end

grad_u=gradient(u,[x,y])

D= 1+x^2+y^2;
for k=1:N 
   D=D+cos(k*pi*x/1).*cos(k*pi*y/1).*param_vec(k);
end

f=-divergence(D*grad_u,[x,y]);