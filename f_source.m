function [f]=f_source(x,y,param_vec_target)
f1=0; f2= x^2 + y^2 + 1; f4=0; f5=-2*x; f7=0; f8=-2*y; N=length(param_vec_target);
for k=1:N
  f1=f1+2*(k*pi)^2*param_vec_target(k)*sin(k*pi*x)*sin(k*pi*y);
  f2=f2+param_vec_target(k)*cos(k*pi*x)*cos(k*pi*y)/N;
  f4=f4+k*param_vec_target(k)*pi*cos(k*pi*x)*sin(k*pi*y);
  f5=f5+k*param_vec_target(k)*pi*cos(k*pi*y)*sin(k*pi*x)/N;
  f7=f7+k*param_vec_target(k)*pi*cos(k*pi*y)*sin(k*pi*x);
  f8=f8+k*param_vec_target(k)*pi*cos(k*pi*x)*sin(k*pi*y)/N;
end
f12=f1*f2;
f45=f4*f5;
f78=f7*f8;
f=f12+f45+f78;


