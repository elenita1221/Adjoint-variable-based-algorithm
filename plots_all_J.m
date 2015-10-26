clear all; % Plots 9 figures
vars={'J3','J4','J5','u_target','Y3','Y4','Y5','Y_target','nruns','N','Mx','My','x','y',...
    'k3','k4','k5','sol3','sol4','sol5','tol','E_u4','E_u4_square','E_u5','E_u5_square',...
    'E_u_target','E_u_target_square','quant3','quant4','quant5'};
load('All_J_N5_50runs-4tol_40x40.mat',vars{:});
E_Y3=mean(Y3(k3,:,:),3)
E_Y4=mean(Y4(k4,:,:),3)
E_Y5=mean(Y5(k5,:,:),3)
E_Y_target=mean(Y_target(:,:),2)'
E_Y3_square=mean(Y3(k3,:,:).*Y3(k3,:,:),3);
E_Y4_square=mean(Y4(k4,:,:).*Y4(k4,:,:),3);
E_Y5_square=mean(Y5(k5,:,:).*Y5(k5,:,:),3);
Var_Y3=E_Y3_square-E_Y3.^2
Var_Y4=E_Y4_square-E_Y4.^2
Var_Y5=E_Y5_square-E_Y5.^2
Var_Y_target=mean(Y_target(:,:).*Y_target(:,:),2)'-(E_Y_target).^2

close all;
if ispc, ms = '\'; else ms = '/'; end
ad = cd;
cd(['..',ms,'Stochastic']);

figure;
max_k=max([k3,k4,k5]); % maximum no. of iterations for plotting on the x-axis
max_quant3=max(quant3(:));max_quant4=max(quant4(:));max_quant5=max(quant5(:));
max_quant=max([max_quant3,max_quant4,max_quant5]);
xlim([1 max_k]); ylim([0 max_quant]);
plot(1:k3,quant3(1:k3),'--g','LineWidth',2);hold on;
plot(1:k4,quant4(1:k4),'-.b','LineWidth',2);hold on;
plot(1:k5,quant5(1:k5),':m','LineWidth',2);hold off;
title(['Mean Convergence in L^2 norm of the estimated solution'])
xlabel('Iterations');ylabel('Quantity of interest ');
legend('||E(u_3)-E(u_{target})||^2','||E(u_4)-E(u_{target})||^2','||E(u_5)-E(u_{target})||^2');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig0_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
max_J3=max(J3(:));max_J4=max(J4(:));max_J5=max(J5(:));max_J=max([max_J3,max_J4,max_J5]);
xlim([1 max_k]); ylim([0 max_J]);
plot(1:k3,J3(1:k3),'--g','LineWidth',2);hold on;
plot(1:k4,J4(1:k4),'-.b','LineWidth',2);hold on;
plot(1:k5,J5(1:k5),':m','LineWidth',2);hold off;
title(['Convergence of the Cost Functional'])
xlabel('Iterations');ylabel('Cost Functional');
legend('J_3','J_4','J_5');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig1_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xlim([1 max_k]);
plot(1:k3,log10(J3(1:k3)),'--g','LineWidth',2);hold on;
plot(1:k4,log10(J4(1:k4)),'-.b','LineWidth',2);hold on;
plot(1:k5,log10(J5(1:k5)),':m','LineWidth',2);hold off;
title(['Convergence of the Log10(Cost functional)'])
xlabel('Iterations');ylabel('Log10(Cost Functional)');
legend('Log10(J_3)','Log10(J_4)','Log10(J_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig2_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for run_counter=1:nruns
 for i=1:Mx+2
  for j=1:My+2
    D_target(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y_target(:,run_counter));
    f_target(i,j,run_counter)=feval('f_source',x(i),y(j),Y_target(:,run_counter));
    D_Y3(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y3(k3,:,run_counter));
    D_Y4(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y4(k4,:,run_counter));
    D_Y5(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y5(k5,:,run_counter));
  end
 end
end
E_D_target=mean(D_target(:,:,:),3);
figure;
surf(x,y,E_D_target);title('Mean of the target diffusion coeff D over the runs')
for i=1:Mx+2
  for j=1:My+2
    D_E_Y3(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y3);  
    D_E_Y4(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y4);
    D_E_Y5(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y5);
  end
end  
figure;
surf(x,y,D_E_Y3(:,:));title('The diffusion coeff D evaluted at the mean of estimated Ys')
%figure;
%surf(x,y,mean(f_target(:,:,:),3));title('Mean of the forcing function f over the runs')
figure;
for run_counter=1:nruns
 plot(x,u_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(x,E_u_target(:,My/2+1),'-r','LineWidth',2);
title('Crossections of the target solution and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig3_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for run_counter=1:nruns
 plot(x,D_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(x,E_D_target(:,My/2+1),'-r','LineWidth',2);
title('Crossections of the target diffusion coeff D and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig4_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(x,E_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(x,D_E_Y3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(x,D_E_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(x,D_E_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Mean of target diffusion coeff versus Mean of estimated diffusion coeff')
legend('E(D_{target})','E(D_3)','E(D_4)','E(D_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig5_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_D_target_square=mean(D_target(:,:,:).*D_target(:,:,:),3);
Var_D_target=E_D_target_square-(E_D_target).^2;
E_D_Y3_square=mean(D_Y3(:,:,:).*D_Y3(:,:,:),3);
Var_D_Y3=E_D_Y3_square-(D_E_Y3).^2;
E_D_Y4_square=mean(D_Y4(:,:,:).*D_Y4(:,:,:),3);
Var_D_Y4=E_D_Y4_square-(D_E_Y4).^2;
E_D_Y5_square=mean(D_Y5(:,:,:).*D_Y5(:,:,:),3);
Var_D_Y5=E_D_Y5_square-(D_E_Y5).^2;
figure;
plot(x,Var_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(x,Var_D_Y3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(x,Var_D_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(x,Var_D_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Variance of target diffusion coeff versus Variance of estimated diffusion coeff')
legend('Var(D_{target})','Var(D_3)','Var(D_4)','Var(D_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig6_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol3_matrix=zeros(Mx+2,My+2,nruns);sol4_matrix=zeros(Mx+2,My+2,nruns);sol5_matrix=zeros(Mx+2,My+2,nruns);
for run_counter=1:nruns
  sol3_matrix(:,:,run_counter)=sol3{run_counter};
  sol4_matrix(:,:,run_counter)=sol4{run_counter};
  sol5_matrix(:,:,run_counter)=sol5{run_counter};
end
E_sol3_matrix=mean(sol3_matrix(:,:,:),3);
E_sol3_square=mean(sol3_matrix(:,:,:).*sol3_matrix(:,:,:),3);
E_sol4_matrix=mean(sol4_matrix(:,:,:),3);%already calculated as E_u4
E_sol5_matrix=mean(sol5_matrix(:,:,:),3);%already calculated as E_u5
 figure;
 surf(x,y,E_sol3_matrix(:,:));title('Mean of the estimated solution over the runs')
 figure;
 surf(x,y,E_u_target(:,:));title('Mean of the target solution over the runs')
figure;
plot(x,E_u_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(x,E_sol3_matrix(:,My/2+1),'--g','LineWidth',2);hold on;
plot(x,E_sol4_matrix(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(x,E_sol5_matrix(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Mean of target solution versus Mean of estimated solution ')
legend('E(u_{target})','E(u_3)','E(u_4)','E(u_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig7_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Var_u_target=E_u_target_square-(E_u_target).^2;  %Variance of u_target
Var_u3=E_sol3_square-(E_sol3_matrix).^2;  %Variance of solution u3
Var_u4=E_u4_square-(E_u4).^2;  %Variance of solution u4
Var_u5=E_u5_square-(E_u5).^2;  %Variance of solution u5
figure;
plot(x,Var_u_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(x,Var_u3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(x,Var_u4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(x,Var_u5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Variance of target solution versus Variance of estimated solution ')
legend('Var(u_{target})','Var(u_3)','Var(u_4)','Var(u_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig8_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_f_target=mean(f_target(:,:,:),3);
figure;
for run_counter=1:nruns
 plot(x,f_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(x,E_f_target(:,My/2+1),'-r','LineWidth',2);
title('Crossections of the target forcing function versus its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig9_',num2str(Mx)]],'-depsc')