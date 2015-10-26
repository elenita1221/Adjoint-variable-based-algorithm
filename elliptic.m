%Omega=[-1, 1, -1, 1];
%Omega=[-pi, pi, -pi, pi];
Omega=[-2,2,-2,2];
Mx=40; My=40;
u_ext='u_exact';
f_src='f_source';
scheme='elliptic5pt';
%calculate the numerical solution at tn=T
tic;
[x,y,dx,dy,u_approx]=feval(scheme,Omega,Mx,My,u_ext,f_src);
toc
%calculate the exact solution at tn=T
u_exact=zeros(Mx+2,My+2);
for i=1:(Mx+2)
for j=1:(My+2)
   u_exact(i,j)=feval(u_ext,x(i),y(j));
end
end
%L^infty error
error_infty=max(max(abs(u_exact-u_approx)))

%L^2 error
error_2= sqrt(dx*dy*sum(sum((u_exact-u_approx).^2)))

%plot the numerical solution
figure;
%[C h]=contour(xcs,ycs,u_approx,10);
%clabel(C,h);
surf(x,y,u_approx)