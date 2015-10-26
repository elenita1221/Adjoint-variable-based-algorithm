clear % Adjoint Elliptic Code using J_3 Cost Functional
N=5; RelError(1)=1e1; nruns=10; eps=1; beta=1e-6; tol=1e-4; eps=2*eps/3; 
Omega=[0, 1, 0, 1]; Mx=40; My=40; %grid size
u_ext='u_exact'; f_src='f_source'; scheme='elliptic5pt';
c=-1;%constant used for plotting the correct number of J functional
mySeed = 10; rng(mySeed); 
 
xes=linspace(Omega(1),Omega(2),Mx+1);
yes=linspace(Omega(3),Omega(4),My+1);
xcs=(xes(1:end-1)+xes(2:end))./2;
ycs=(yes(1:end-1)+yes(2:end))./2;
x=[xes(1),xcs,xes(Mx+1)];y=[yes(1),ycs,yes(Mx+1)];
for run_counter=1:nruns 
     %Y_target(1:N,run_counter) = 0.4 + (0.6-0.4).*rand(1,N); 
     Y_target(1:N,run_counter) = rand(1,N);
   for m=1:(Mx+2)
   for n=1:(My+2)
     u_target(m,n,run_counter)=feval(u_ext,x(m),y(n),N,Y_target(:,run_counter));
   end
   end
end

for i=1:Mx
 for j=1:My
    cos_prod(i,j,1:N)=cos((1:N)*pi*x(i+1)).*cos((1:N)*pi*y(j+1))./N;
    cos_sum(i,j)=sum(cos_prod(i,j,:))/N;
 end
end  

k=1; sol=cell(1,nruns);
for run_counter=1:nruns 
%Y_ini=0.4 + (0.6-0.4).*rand(1,N); %initial guess for 1st iteration
Y_ini = rand(1,N); 
Y(k,1:N,run_counter)=Y_ini;
[x,y,dx,dy,u,int_k,params]=feval(scheme,N,Y(k,1:N,run_counter),Y_target(1:N,run_counter),Omega,Mx,My,u_ext,f_src);%SE elliptic solver
sol{run_counter}=u;
J(run_counter,k) =sum(sum(((u(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2; %cost functional
end
J_new(k)=mean(J(1:nruns,k));

while (RelError(k)>tol) 
     eps=3*eps/2; 
     k=k+1;
 for run_counter=1:nruns 
     u=sol{run_counter};
     [~,~,~,~,eta,~,~]=feval('elliptic5pt_adjoint',N,Y(k-1,1:N,run_counter),Omega,Mx,My,u-u_target(:,:,run_counter));%AE elliptic solver
     gradx_u = [diff(u(:,:),1,2)./dx];
     grady_u = [diff(u(:,:),1,1)./dy];
     gradx_eta=[diff(eta(:,:),1,2)./dx];
     grady_eta=[diff(eta(:,:),1,1)./dy];
     matrix1=(gradx_u(:,1:Mx).*gradx_eta(:,1:Mx)+gradx_u(:,2:Mx+1).*gradx_eta(:,2:Mx+1))/2;
     matrix2=(grady_u(1:My,:).*grady_eta(1:My,:)+grady_u(2:My+1,:).*grady_eta(2:My+1,:))/2;
     temp1=-sum(sum(cos_sum(:,:).*matrix1(2:My+1,:)))*dx*dy;
     temp2=-sum(sum(cos_sum(:,:).*matrix2(:,2:Mx+1)))*dx*dy;
     for i=1:N
       temp3=beta*sum(sum(params.D(:,:).*cos_prod(:,:,i)))*dx*dy;
       dJdY(k-1,i,run_counter)=temp1+temp2+temp3;
     end
  Y(k,1:N,run_counter)=Y(k-1,1:N,run_counter)-eps*dJdY(k-1,1:N,run_counter);
  [x,y,dx,dy,u,int_k,params]=feval(scheme,N,Y(k,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
  sol{run_counter}=u;
  J(run_counter,k)=sum(sum(((u(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2;
end
J_new(k)=mean(J(1:nruns,k));
 
  while J_new(k)>=J_new(k-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end; 
    for run_counter=1:nruns
     Y(k,1:N,run_counter)=Y(k-1,1:N,run_counter)-eps*dJdY(k-1,1:N,run_counter);
     [x,y,dx,dy,u,int_k,params]=feval(scheme,N,Y(k,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol{run_counter}=u;
     J(run_counter,k)=sum(sum(((u(1:Mx+2,1:My+2)-u_target(1:Mx+2,1:My+2,run_counter)).^2).*dx.*dy))/2+beta*int_k/2;
    end 
    J_new(k)=mean(J(1:nruns,k));
  end
  
RelError(k)=abs(J_new(k)-J_new(k-1))/abs(J_new(k));
end


filename=strcat('J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol_',num2str(Mx),'x',num2str(My),'.mat')
save(filename) 

    
    
    
    
    
