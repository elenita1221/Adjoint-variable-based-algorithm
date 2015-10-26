function  [x,y,dx,dy,u,int_D,params]= elliptic5pt_adjoint(N,param_vec,Omega,Mx,My,f_src)
%	Dirichlet boundary conditions
%Input:
%  Omega      -- [x1,x2]x[y1,y2] computational domain
%  Mx,My      -- number of cells in x,y (M+1 grid points)
%  u_ext      -- function, exact solution u(x,y), to get I.C. 
%Output:
%  xes,yes    -- (x(i),y(j)) coordinates of grid point (i,j)
%  u          -- grid funtion at t=T	
%Local:
%  dx,dy      -- grid size in x, y directions
%  A          -- global matrix
%  vec,rhs    -- global solution vector and rhs vector

xes=linspace(Omega(1),Omega(2),Mx+1);
yes=linspace(Omega(3),Omega(4),My+1);
xcs=(xes(1:end-1)+xes(2:end))./2;
ycs=(yes(1:end-1)+yes(2:end))./2;
x=[xes(1),xcs,xes(Mx+1)];y=[yes(1),ycs,yes(Mx+1)];
dx=(Omega(2)-Omega(1))./Mx;
dy=(Omega(4)-Omega(3))./My;
u=zeros(Mx+2,My+2);
vec=zeros(Mx*My,1);
rhs=zeros(Mx*My,1);
A=sparse(Mx*My,Mx*My);
cx=1/dx^2;
cy=1/dy^2;

%boundary condition
for i=1:Mx+2
for j=1:My+2
  if(glbidx(i,j,Mx,My)==0)  %boundary points
    u(i,j)=0;
  end
end
end

for i=1:Mx+1
 for j=1:My
    params.Dx(i,j)=feval('Dcoeff',xes(i),ycs(j),N,param_vec);
 end
end

for i=1:Mx
 for j=1:My+1
    params.Dy(i,j)=feval('Dcoeff',xcs(i),yes(j),N,param_vec);
 end
end 

for i=1:Mx
 for j=1:My
    params.D(i,j)=feval('Dcoeff',xcs(i),ycs(j),N,param_vec);
 end
end 


int_D=sum(params.D(:).*params.D(:).*dx.*dy);
%Assembling the global matrix and RHS vector-------------------------------
for i=1:Mx+2
for j=1:My+2
  m1=glbidx(i,j,Mx,My);
  if(m1==0)
    continue;
  end
  rhs(m1)=f_src(i,j);%RHS vector
  rhs(m1)=-rhs(m1);
  
% Global matrix -----------------------------------------------------------
 if (i==2) && (j==2)
   %diagonal entry (i,j)
   A(m1,m1)=-(params.Dx(i,j-1)*cx+2*params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
   +2*params.Dy(i-1,j-1)*cy);
     i2=i+1;%(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     i2=i; %(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
     rhs(m1)=rhs(m1)-cx*u(i-1,j)*2*params.Dx(i-1,j-1)-cy*u(i,j-1)*2*params.Dy(i-1,j-1);
     continue;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if (i>2) && (i<=Mx) && (j==2)
   %diagonal entry (i,j)
   A(m1,m1)=-(params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
   +2*params.Dy(i-1,j-1)*cy);
     i2=i-1;%(i-1,j)  
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     i2=i+1;%(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     i2=i;%(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
     rhs(m1)=rhs(m1)-cy*u(i,j-1)*2*params.Dy(i-1,j-1);
     continue;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==(Mx+1)) && (j==2)
   %diagonal entry (i,j)
   A(m1,m1)=-(2*params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
   +2*params.Dy(i-1,j-1)*cy);
     i2=i-1;%(i-1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     i2=i;%(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
     rhs(m1)=rhs(m1)-cx*u(i+1,j)*2*params.Dx(i,j-1)-cy*u(i,j-1)*2*params.Dy(i-1,j-1);
     continue;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==(Mx+1)) && (j>2) && (j<=My)
   %diagonal entry (i,j)
   A(m1,m1)=-(2*params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
   +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i-1;%(i-1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     i2=i;%(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
     rhs(m1)=rhs(m1)-cx*u(i+1,j)*2*params.Dx(i,j-1);
     continue;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==(Mx+1)) && (j==(My+1))
   %diagonal entry (i,j)
   A(m1,m1)=-(2*params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+2*params.Dy(i-1,j)*cy...
   +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i-1;%(i-1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     rhs(m1)=rhs(m1)-cx*u(i+1,j)*2*params.Dx(i,j-1)-cy*u(i,j+1)*2*params.Dy(i-1,j);
     continue;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (i>2) && (i<=Mx) && (j==(My+1))
   %diagonal entry (i,j)
   A(m1,m1)=-(params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+2*params.Dy(i-1,j)*cy...
   +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i-1;%(i-1,j)  
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     i2=i+1;%(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     rhs(m1)=rhs(m1)-cy*u(i,j+1)*2*params.Dy(i-1,j);
     continue;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==2) && (j==(My+1))
   %diagonal entry (i,j)
   A(m1,m1)=-(params.Dx(i,j-1)*cx+2*params.Dx(i-1,j-1)*cx+2*params.Dy(i-1,j)*cy...
   +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i+1; %(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     rhs(m1)=rhs(m1)-cx*u(i-1,j)*2*params.Dx(i-1,j-1)-cy*u(i,j+1)*2*params.Dy(i-1,j);
     continue;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==2) && (j>2) && (j<=My)
   %diagonal entry (i,j)
   A(m1,m1)=-(params.Dx(i,j-1)*cx+2*params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
   +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i+1; %(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     i2=i;%(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
     rhs(m1)=rhs(m1)-cx*u(i-1,j)*2*params.Dx(i-1,j-1);
     continue;
end
%Interior Points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A(m1,m1)=-(params.Dx(i,j-1)*cx+params.Dx(i-1,j-1)*cx+params.Dy(i-1,j)*cy...
 +params.Dy(i-1,j-1)*cy);
     i2=i;%(i,j-1)
     j2=j-1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j-1)*cy;
     i2=i-1;%(i-1,j)  
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i-1,j-1)*cx;
     i2=i+1; %(i+1,j)
     j2=j;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dx(i,j-1)*cx;
     i2=i;%(i,j+1)
     j2=j+1;   	 	
     m2=glbidx(i2,j2,Mx,My);
     A(m1,m2)=params.Dy(i-1,j)*cy;
end
end
%eig(A)
%vec=(A)\(rhs);
A=-A;	%make sure the matrix sent to CG solver is SPD
rhs=-rhs;
[vec,flag]=pcg(A,rhs,1.e-6,500);

for i=2:(Mx+1)
for j=2:(My+1)
  m=glbidx(i,j,Mx,My);
  u(i,j)=vec(m);
end
end


function m=glbidx(i,j,Mx,My)
% subfuction to calculate the index of u(i,j) in the gobal vector
if(i==1||j==1||i==Mx+2||j==My+2)
  m=0;
else
  m=(i-1)+(j-2)*Mx;
end