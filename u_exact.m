function u = u_exact(x,y,N,param_vec)
%U_EXACT exact solution 
u=0;
for k=1:N 
   u=u+sin(k*pi*x/1).*sin(k*pi*y/1).*param_vec(k);
end
