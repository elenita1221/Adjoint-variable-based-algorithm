clear all %Adjoint Elliptic Code using J3 or J4 OR J5 Cost Functional
N=5; RelError(1)=1e1; nruns=50; eps=1; beta=1e-6; tol=1e-4; eps=2*eps/3; 
Omega=[0, 1, 0, 1]; Mx=40; My=40; %grid size
u_ext='u_exact'; f_src='f_source'; scheme='elliptic5pt'; 
E_u_target=zeros(Mx+2,My+2);E_k_square=zeros(Mx,My);
E_u_target_square=zeros(Mx+2,My+2);
mySeed = 10; rng(mySeed); 

xes=linspace(Omega(1),Omega(2),Mx+1);
yes=linspace(Omega(3),Omega(4),My+1);
xcs=(xes(1:end-1)+xes(2:end))./2;
ycs=(yes(1:end-1)+yes(2:end))./2;
x=[xes(1),xcs,xes(Mx+1)];y=[yes(1),ycs,yes(Mx+1)];
for run_counter=1:nruns 
     %Y_target(1:N,run_counter) = 0.4 + (0.6-0.4).*rand(1,N); 
     Y_target(1:N,run_counter) = rand(1,N); 
     %Y_target(1:N,run_counter) = 0.5+(1:N)/100+run_counter/1000;
   for m=1:(Mx+2)
   for n=1:(My+2)
     u_target(m,n,run_counter)=feval(u_ext,x(m),y(n),N,Y_target(:,run_counter));
   end
   end
   E_u_target=E_u_target+u_target(:,:,run_counter);
   E_u_target_square=E_u_target_square+u_target(:,:,run_counter).*u_target(:,:,run_counter);
end
E_u_target=E_u_target/nruns; E_u_target_square=E_u_target_square/nruns;

for i=1:Mx
 for j=1:My
    cos_prod(i,j,1:N)=cos((1:N)*pi*x(i+1)).*cos((1:N)*pi*y(j+1))./N;
    cos_sum(i,j)=sum(cos_prod(i,j,:))/N;
 end
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k3=1; sol3=cell(1,nruns);  %J3 cost functional
E_u3=zeros(Mx+2,My+2); E_u3_square=zeros(Mx+2,My+2);
for run_counter=1:nruns 
%Y_ini=0.4 + (0.6-0.4).*rand(1,N); %initial guess for 1st iteration
Y_ini(1:N,run_counter) = rand(1,N);
Y3(k3,1:N,run_counter)=Y_ini(1:N,run_counter);
[x,y,dx,dy,u3,int_k,params]=feval(scheme,N,Y3(k3,1:N,run_counter),Y_target(1:N,run_counter),Omega,Mx,My,u_ext,f_src);%SE elliptic solver
sol3{run_counter}=u3;
J(run_counter,k3) =sum(sum(((u3(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2; %cost functional
E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3;
end
J3(k3)=mean(J(1:nruns,k3));
E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
quant3(k3)=sum(sum(((E_u3(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));

while (RelError(k3)>tol) 
     eps=3*eps/2; 
     k3=k3+1;
     E_u3=zeros(Mx+2,My+2); E_u3_square=zeros(Mx+2,My+2);
 for run_counter=1:nruns 
     u3=sol3{run_counter};
     [~,~,~,~,eta3,~,~]=feval('elliptic5pt_adjoint',N,Y3(k3-1,1:N,run_counter),Omega,Mx,My,u3-u_target(:,:,run_counter));%AE elliptic solver
     gradx_u = [diff(u3(:,:),1,2)./dx];
     grady_u = [diff(u3(:,:),1,1)./dy];
     gradx_eta=[diff(eta3(:,:),1,2)./dx];
     grady_eta=[diff(eta3(:,:),1,1)./dy];
     matrix1=(gradx_u(:,1:Mx).*gradx_eta(:,1:Mx)+gradx_u(:,2:Mx+1).*gradx_eta(:,2:Mx+1))/2;
     matrix2=(grady_u(1:My,:).*grady_eta(1:My,:)+grady_u(2:My+1,:).*grady_eta(2:My+1,:))/2;
     temp1=-sum(sum(cos_sum(:,:).*matrix1(2:My+1,:)))*dx*dy;
     temp2=-sum(sum(cos_sum(:,:).*matrix2(:,2:Mx+1)))*dx*dy;
     for i=1:N
       temp3=beta*sum(sum(params.D(:,:).*cos_prod(:,:,i)))*dx*dy;
       dJdY(k3-1,i,run_counter)=temp1+temp2+temp3;
     end
  Y3(k3,1:N,run_counter)=Y3(k3-1,1:N,run_counter)-eps*dJdY(k3-1,1:N,run_counter);
  [x,y,dx,dy,u3,int_k,params]=feval(scheme,N,Y3(k3,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
  sol3{run_counter}=u3;
  J(run_counter,k3)=sum(sum(((u3(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2;
  E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3; 
end
J3(k3)=mean(J(1:nruns,k3));
E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
quant3(k3)=sum(sum(((E_u3(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));

  while J3(k3)>=J3(k3-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u3=zeros(Mx+2,My+2); 
    E_u3_square=zeros(Mx+2,My+2);
    for run_counter=1:nruns
     Y3(k3,1:N,run_counter)=Y3(k3-1,1:N,run_counter)-eps*dJdY(k3-1,1:N,run_counter);
     [x,y,dx,dy,u3,int_k,params]=feval(scheme,N,Y3(k3,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol3{run_counter}=u3;
     J(run_counter,k3)=sum(sum(((u3(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2;
     E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3; 
    end 
    J3(k3)=mean(J(1:nruns,k3));
    E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
    quant3(k3)=sum(sum(((E_u3(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
  end
  
RelError(k3)=abs(J3(k3)-J3(k3-1))/abs(J3(k3));
end  
disp('J3 finished!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dJdY
eps=1; eps=2*eps/3;
k4=1; sol4=cell(1,nruns); c=0;  %J4 cost functional
E_u4=zeros(Mx+2,My+2); E_u4_square=zeros(Mx+2,My+2);
for run_counter=1:nruns 
%Y_ini=0.4 + (0.6-0.4).*rand(1,N); %initial guess for 1st iteration
Y4(k4,1:N,run_counter)=Y_ini(1:N,run_counter);
[x,y,dx,dy,u4,int_k,params]=feval(scheme,N,Y4(k4,1:N,run_counter),Y_target(1:N,run_counter),Omega,Mx,My,u_ext,f_src);%SE elliptic solver
sol4{run_counter}=u4;
E_u4=E_u4+u4; E_u4_square=E_u4_square+u4.*u4; E_k_square=E_k_square+params.D.*params.D;
end;
E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
E_u4_square=E_u4_square/nruns;
quant4(k4)=sum(sum(((E_u4(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
J4(k4)=quant4(k4)/2+...
    +c*sum(sum(((E_u4_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
    +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional

old_E_u4=E_u4; old_E_u4_square=E_u4_square;
while RelError(k4)>tol
     eps=3*eps/2;
     k4=k4+1;
     E_u4=zeros(Mx+2,My+2); E_k_square=zeros(Mx,My);
     E_u4_square=zeros(Mx+2,My+2);
 for run_counter=1:nruns 
    u4=sol4{run_counter};
    [~,~,~,~,eta4,~,~]=feval('elliptic5pt_adjoint',N,Y4(k4-1,1:N,run_counter),Omega,Mx,My,old_E_u4-E_u_target+c*u4.*(old_E_u4_square-E_u_target_square));%AE elliptic solver
    gradx_u = [diff(u4(:,:),1,2)./dx];
    grady_u = [diff(u4(:,:),1,1)./dy];
    gradx_eta=[diff(eta4(:,:),1,2)./dx];
    grady_eta=[diff(eta4(:,:),1,1)./dy];
    matrix1=(gradx_u(:,1:Mx).*gradx_eta(:,1:Mx)+gradx_u(:,2:Mx+1).*gradx_eta(:,2:Mx+1))/2;
    matrix2=(grady_u(1:My,:).*grady_eta(1:My,:)+grady_u(2:My+1,:).*grady_eta(2:My+1,:))/2;
    temp1=-sum(sum(cos_sum(:,:).*matrix1(2:My+1,:)))*dx*dy;
    temp2=-sum(sum(cos_sum(:,:).*matrix2(:,2:Mx+1)))*dx*dy;
    for i=1:N
       temp3=beta*sum(sum(params.D(:,:).*cos_prod(:,:,i)))*dx*dy;
       dJdY(k4-1,i,run_counter)=temp1+temp2+temp3;
    end
    Y4(k4,1:N,run_counter)=Y4(k4-1,1:N,run_counter)-eps*dJdY(k4-1,1:N,run_counter);
    [x,y,dx,dy,u4,int_k,params]=feval(scheme,N,Y4(k4,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
    sol4{run_counter}=u4;
    E_u4=E_u4+u4; E_k_square=E_k_square+params.D.*params.D;
    E_u4_square=E_u4_square+u4.*u4; 
 end
 E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
 E_u4_square=E_u4_square/nruns;
 quant4(k4)=sum(sum(((E_u4(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
 J4(k4)=quant4(k4)/2+...
    +c*sum(sum(((E_u4_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
    +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional

  while J4(k4)>=J4(k4-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u4=zeros(Mx+2,My+2); E_k_square=zeros(Mx,My);
    E_u4_square=zeros(Mx+2,My+2);
    for run_counter=1:nruns
     Y4(k4,1:N,run_counter)=Y4(k4-1,1:N,run_counter)-eps*dJdY(k4-1,1:N,run_counter);
     [x,y,dx,dy,u4,int_k,params]=feval(scheme,N,Y4(k4,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol4{run_counter}=u4;
     E_u4=E_u4+u4; E_k_square=E_k_square+params.D.*params.D;
     E_u4_square=E_u4_square+u4.*u4; 
    end;
    E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
    E_u4_square=E_u4_square/nruns;
    quant4(k4)=sum(sum(((E_u4(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
    J4(k4)=quant4(k4)/2+...
     +c*sum(sum(((E_u4_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
     +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional
  end
old_E_u4=E_u4; old_E_u4_square=E_u4_square;
RelError(k4)=abs(J4(k4)-J4(k4-1))/abs(J4(k4));
end
disp('J4 finished!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dJdY E_k_square
eps=1; eps=2*eps/3;
k5=1; sol5=cell(1,nruns); c=1;  %J5 cost functional
E_u5=zeros(Mx+2,My+2); E_u5_square=zeros(Mx+2,My+2); E_k_square=zeros(Mx,My);
for run_counter=1:nruns 
%Y_ini=0.4 + (0.6-0.4).*rand(1,N); %initial guess for 1st iteration 
Y5(k5,1:N,run_counter)=Y_ini(1:N,run_counter);
[x,y,dx,dy,u5,int_k,params]=feval(scheme,N,Y5(k5,1:N,run_counter),Y_target(1:N,run_counter),Omega,Mx,My,u_ext,f_src);%SE elliptic solver
sol5{run_counter}=u5;
E_u5=E_u5+u5; E_u5_square=E_u5_square+u5.*u5; E_k_square=E_k_square+params.D.*params.D;
end;
E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
E_u5_square=E_u5_square/nruns;
quant5(k5)=sum(sum(((E_u5(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
J5(k5)=quant5(k5)/2+...
    +c*sum(sum(((E_u5_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
    +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional

old_E_u5=E_u5; old_E_u5_square=E_u5_square;
while RelError(k5)>tol
     eps=3*eps/2;
     k5=k5+1;
     E_u5=zeros(Mx+2,My+2); E_k_square=zeros(Mx,My);
     E_u5_square=zeros(Mx+2,My+2);
 for run_counter=1:nruns 
    u5=sol5{run_counter};
    [~,~,~,~,eta5,~,~]=feval('elliptic5pt_adjoint',N,Y5(k5-1,1:N,run_counter),Omega,Mx,My,old_E_u5-E_u_target+c*u5.*(old_E_u5_square-E_u_target_square));%AE elliptic solver
    gradx_u = [diff(u5(:,:),1,2)./dx];
    grady_u = [diff(u5(:,:),1,1)./dy];
    gradx_eta=[diff(eta5(:,:),1,2)./dx];
    grady_eta=[diff(eta5(:,:),1,1)./dy];
    matrix1=(gradx_u(:,1:Mx).*gradx_eta(:,1:Mx)+gradx_u(:,2:Mx+1).*gradx_eta(:,2:Mx+1))/2;
    matrix2=(grady_u(1:My,:).*grady_eta(1:My,:)+grady_u(2:My+1,:).*grady_eta(2:My+1,:))/2;
    temp1=-sum(sum(cos_sum(:,:).*matrix1(2:My+1,:)))*dx*dy;
    temp2=-sum(sum(cos_sum(:,:).*matrix2(:,2:Mx+1)))*dx*dy;
    for i=1:N
       temp3=beta*sum(sum(params.D(:,:).*cos_prod(:,:,i)))*dx*dy;
       dJdY(k5-1,i,run_counter)=temp1+temp2+temp3;
    end
    Y5(k5,1:N,run_counter)=Y5(k5-1,1:N,run_counter)-eps*dJdY(k5-1,1:N,run_counter);
    [x,y,dx,dy,u5,int_k,params]=feval(scheme,N,Y5(k5,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
    sol5{run_counter}=u5;
    E_u5=E_u5+u5; E_k_square=E_k_square+params.D.*params.D;
    E_u5_square=E_u5_square+u5.*u5; 
 end
 E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
 E_u5_square=E_u5_square/nruns;
 quant5(k5)=sum(sum(((E_u5(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
 J5(k5)=quant5(k5)/2+...
    +c*sum(sum(((E_u5_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
    +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional

  while J5(k5)>=J5(k5-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u5=zeros(Mx+2,My+2); E_k_square=zeros(Mx,My);
    E_u5_square=zeros(Mx+2,My+2);
    for run_counter=1:nruns
     Y5(k5,1:N,run_counter)=Y5(k5-1,1:N,run_counter)-eps*dJdY(k5-1,1:N,run_counter);
     [x,y,dx,dy,u5,int_k,params]=feval(scheme,N,Y5(k5,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol5{run_counter}=u5;
     E_u5=E_u5+u5; E_k_square=E_k_square+params.D.*params.D;
     E_u5_square=E_u5_square+u5.*u5; 
    end;
    E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
    E_u5_square=E_u5_square/nruns;
    quant5(k5)=sum(sum(((E_u5(1:Mx+2,1:My+2)-E_u_target(1:Mx+2,1:My+2)).^2).*dx.*dy));
    J5(k5)=quant5(k5)/2+...
     +c*sum(sum(((E_u5_square(1:Mx+2,1:My+2)-E_u_target_square(1:Mx+2,1:My+2)).^2).*dx.*dy))/4+...
     +beta*sum(E_k_square(:).*dx.*dy)/2; %cost functional
  end
old_E_u5=E_u5; old_E_u5_square=E_u5_square;
RelError(k5)=abs(J5(k5)-J5(k5-1))/abs(J5(k5));
end

filename=strcat('All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol_',num2str(Mx),'x',num2str(My),'.mat')
save(filename)
