function [ecs,int_D] = ec_model1(N,param_vec,ecfinal,ets,g,f,varargin)


%
%  param_vec-vector of parameters, currently param_vec = [D,kp], just two
%    though it can be adjusted as need be so long as you adjust the program
%    accordingly, in particular the subprogram advance_ecs would need to be
%    changed.
%  ecini-initial epithelial cell density in the grid
%  ets-experimental times at which we want to know the epithelial cell 
%    density in the grid
%  g-grid info including
%    g.nx-number of cells in x-direction
%    g.ny-number of cells in y-direction
%    g.dx-size of cells in x-direction (we assume uniform spacing)
%    g.dy-size of cells in y-direction (again uniform spacing assumed)
%  Note:  The info inside g is later augmented with the appropriate time
%    step to use during each step in the numerical process (changes with
%    time).

%%  Extract info from param_vec
%  D-diffusion coefficient, corresponds to diffusion rate of cells
%  kp-growth or proliferation rate of epithelial cells
%  hc-hill coefficient, power in the hill term in the diffusion
%  Note:  To make the hill coefficient a free parameter, you need to swap
%  the commented and uncommented lines below for hc.  You would also need
%  to make the right changes in the MCMC program that called this program.
% params.D = param_vec;
[xM,yM] = meshgrid(g.xcs,g.ycs);
params.D= ones(g.ny,g.nx)+xM.^2+yM.^2;
for i=1:N 
     params.D=params.D+cos(i*pi*xM/0.1).*cos(i*pi*yM/0.07).*param_vec(i);
end
params.D=params.D*(1e-3);
 % params.kp = param_vec(2);
  % params.hc = param_vec(3);
 % params.hc = 2;

  %%  Consts for good numerical solving and other options
  ecm.safety_net =.5;
  ecm.plot_soln = 1;
  beta=1e-3; 
int_D=beta*sum(params.D(:).^2.*g.dx.*g.dy)/2; 
  
  for vac = 1:2:length(varargin)
    ecm.(varargin{vac}) = varargin{vac+1};
  end

  %%  Get maximum dt allowable for stability
  %  This includes a safety buffer fraction of 0.9.
  temp_quant = g.dx^2*g.dy^2/(max(params.D(:))*(g.dx^2+g.dy^2));
  maxdt = ecm.safety_net*temp_quant;
  
  ecs = zeros([size(ecfinal),length(ets)]);
  ecs(:,:,end) = ecfinal;
  my_plot(ecm.plot_soln,ecfinal)

  for etc = length(ets)-1:-1:1

    exp_dt = ets(etc+1)-ets(etc);
    nsteps = ceil(exp_dt/maxdt);
    comp_dt = exp_dt/nsteps;
    g.dt = abs(comp_dt);

    ecs(:,:,etc) = ecs(:,:,etc+1);
    for ctc = 1:nsteps
      ecs(:,:,etc) = advance_ecs(ecs(:,:,etc),g,params,f(:,:,etc));
    end

    my_plot(ecm.plot_soln,ecs(:,:,etc))

  end
end

function my_plot(plot_soln,ecs)
  if plot_soln
    surf(ecs);
    zlim([0,1]);
    hold on; 
    contour(ecs,[0.5 0.5]);
    hold off;
    pause;
  end
end

function ecs = advance_ecs(ecs,g,params,f)
  
 % temp = ecs.^params.hc;
 % nonlinpartofdiffcoef = temp./(temp+(1-ecs).^params.hc);
  
  grad_x = [(ecs(:,1)-0)./(g.dx/2),...
    diff(ecs,1,2)./g.dx,...
    (0-ecs(:,end))./(g.dx/2)];
%  diff_coef_x = [ones(size(nonlinpartofdiffcoef(1,:)));...
%    (nonlinpartofdiffcoef(1:end-1,:)+nonlinpartofdiffcoef(2:end,:))./2;...
%    ones(size(nonlinpartofdiffcoef(end,:)))];
    
  grad_y = [(ecs(1,:)-0)./(g.dy/2);...
    diff(ecs,1,1)./g.dy;...
    (0-ecs(end,:))./(g.dy/2)];
%  diff_coef_y = [(1+nonlinpartofdiffcoef(:,1))./2,...
%    (nonlinpartofdiffcoef(:,1:end-1)+nonlinpartofdiffcoef(:,2:end))./2,...
%    (1+nonlinpartofdiffcoef(:,end))./2];
  
  ecs = ecs+g.dt*(params.D.*(diff(grad_x,1,2)./g.dx+...
    diff(grad_y,1,1)./g.dy)+f);

end
