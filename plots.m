clear; % Plots 12 figures
vars={'J_new','u_target','Y','Y_target','nruns','N','Mx','My','k','x','y','sol','u','c','tol'};
load('J4_N5_10runs_0.0001tol_40x40.mat',vars{:});
E_Y=mean(Y(k,:,:),3)
E_Y_target=mean(Y_target(:,:),2)'
E_Y_square=mean(Y(k,:,:).*Y(k,:,:),3);
Var_Y=E_Y_square-E_Y.^2
%Var_u=old_E_u_square-(old_E_u).^2;  %Variance of the solution u
u_target_mean=mean(u_target(:,:,:),3);

close all;
if ispc, ms = '\'; else ms = '/'; end
ad = cd;
cd(['..',ms,'Stochastic']);

figure;
plot(1:k,J_new(1:k),'LineWidth',2);
title(['Convergence of the cost functional J_',num2str(c+4)])
xlabel('Iterations');ylabel('Cost Functional');
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig1_',num2str(Mx)]],'-depsc')

figure;
plot(1:k,log10(J_new(1:k)),'LineWidth',2);
title(['Convergence of the Log10 of cost functional J_',num2str(c+4)])
xlabel('Iterations');ylabel('Log10(Cost Functional)');
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig2_',num2str(Mx)]],'-depsc')

for run_counter=1:nruns
 for i=1:Mx+2
  for j=1:My+2
    D(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y_target(:,run_counter));
    f(i,j,run_counter)=feval('f_source',x(i),y(j),Y_target(:,run_counter));
  end
 end
end
D_mean=mean(D(:,:,:),3);
%figure;
%surf(x,y,D_mean);title('Mean of the target diffusion coefficient D over the runs')
for i=1:Mx+2
  for j=1:My+2
    D_E_Y(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y);
  end
end  
%figure;
%surf(x,y,D_E_Y(:,:));title('The diffusion coefficient D evaluted at the mean of estimated Ys')

%figure;
%surf(x,y,mean(f(:,:,:),3));title('Mean of the forcing function f over the runs')
figure;
for run_counter=1:nruns
 plot(x,u_target(:,My/2+1,run_counter),'b');
 hold on;
end
plot(x,u_target_mean(:,My/2+1),'r','LineWidth',2);
title('Crossections of the target solution and its mean')
hold off;
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig3_',num2str(Mx)]],'-depsc')

figure;
for run_counter=1:nruns
 plot(x,D(:,My/2+1,run_counter),'b');
 hold on;
end
plot(x,D_mean(:,My/2+1),'r','LineWidth',2);
title('Crossections of the target diffusion coeff D versus its mean')
hold off;
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig4_',num2str(Mx)]],'-depsc')

figure;
plot(x,D_mean(:,My/2+1),'r','LineWidth',2);hold on;
plot(x,D_E_Y(:,My/2+1),'b','LineWidth',2);
title('Crossections of mean of target diffusion coeff versus mean of estimated diffusion coeff')
hold off;
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig5_',num2str(Mx)]],'-depsc')

figure;
for run_counter=1:nruns
 plot(x,f(:,My/2+1,run_counter),'b');
 hold on;
end
plot(x,mean(f(:,My/2+1,:),3),'r','LineWidth',2);
title('Crossections of the target forcing function versus its mean')
hold off;
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig6_',num2str(Mx)]],'-depsc')

% figure;
 sol_matrix=zeros(Mx+2,My+2,nruns);
 for run_counter=1:nruns
  sol_matrix(:,:,run_counter)=sol{run_counter};
 end
 sol_matrix_mean=mean(sol_matrix(:,:,:),3);
% surf(x,y,sol_matrix_mean(:,:));title('Mean of the estimated solution over the runs')
% figure;
% surf(x,y,u_target_mean(:,:));title('Mean of the target solution over the runs')

figure;
plot(x,u_target_mean(:,My/2+1),'r');
hold on;
plot(x,sol_matrix_mean(:,My/2+1),'b','LineWidth',2);
title('Crossections of mean of target solution versus mean of estimated solution ')
hold off;
print(gcf,[['Figures_J',num2str(c+4),'_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig7_',num2str(Mx)]],'-depsc')

