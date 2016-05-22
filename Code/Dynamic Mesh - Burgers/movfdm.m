function [t,x,u]=movfdm(npde,jmax,tspan,x0,xt0,u0,vt0,PDE_F,PDE_G,BC_L,BC_R, ...
                        TauFun,monitor,reltol,abstol,MonitorFun,mmpde,alpha_type,alpha)
%
% [t,x,u]=movfdm(npde,jmax,tspan,x0,xt0,u0,vt0,PDE_F,PDE_G,BC_L,BC_R, ...
%                TauFun,monitor,reltol,abstol,MonitorFun,mmpde,alpha_type,alpha)
%
% This function integrates a system of (npde) one-dimensional PDEs in the form
%
%   F(t, x, u, u_x, u_t) = (d/dx) G (t, x, u, u_x, u_t)  for a < x < b
%                      0 = B_L(t, x, u, u_x, u_t)         at x = a
%                      0 = B_R(t, x, u, u_x, u_t)         at x = b
%
% using an MMPDE moving mesh strategy. The modified MMPDE5 is used as
% the default strategy.
%  
% The physical and mesh PDEs are discretized using central finite differences and
% the resulting ODE system is integrated using the Matlab built-in solver ode15i.
%
% Note that the first 11 arguments must be specified by the user. The other 8
% arguments (TauFun,monitor,reltol,abstol,MonitorFun,mmpde,alpha_type,alpha) are optional.
%
% output variables:
%
% t:        column vector of time instants, t(1:N)
% x:        meshes at all time levels, x(j,n), j=1:jmax, n=1:N
% u:        solutions at all time levels, u(i,j,n), i=1:npde, j=1:jmax, n=1:N
%
% input variables:
%
% npde:     the number of physical PDEs contained in the system to be integrated.
% jmax:     the number of mesh points.
% tspan:    an array specifying the time interval for integration.
% x0:       column vector of size jmax, initial mesh.
% xt0:      column vector of size jmax, approximate values for initial mesh speed
%           which are required by the implicit DAE solver ode15i.
% u0:       array of size npde by jmax, intial values for the physical solution.
%           u0(i,j): i=1:npde, j=1:jmax
% vt0:      array of size npde by jmax, approximate values for the initial time
%           derivatives of the physical solution which are required by the implicit
%           DAE solver ode15i.
% PDE_F,PDE_G,BC_L,BC_R: functions used to define the system of physical PDEs:
%
%           PDE_F := F(t,x,u,u_x,u_t)
%           PDE_G := G(t,x,u,u_x,u_t)
%            BC_L := B_L(t,x,u,u_x,u_t)
%            BC_R := B_R(t,x,u,u_x,u_t)
%
%    function f=PDE_F(t,x,u,ux,ut)
%    function g=PDE_G(t,x,u,ux,ut)
%
%           input variables:
%           t is scalar, x=x(1:m) is a row vector, and u=u(1:npde,1:m), ux=ux(1:npde,1:m),
%           and ut=ut(1:npde,1:m) are matrices.
%
%           output variable:
%           f=f(1:npde,1:m) and g=g(1:npde,1:m) should be matrices.
%
%    function res=BC_L(t,x,u,ux,ut)
%    function res=BC_R(t,x,u,ux,ut)
%
%           input variables:
%           t and x are scalar and u=u(1:npde), ux=ux(1:npde), and ut=ut(1:npde) are
%           column vectors.
%
%           output variable:
%           res=res(1:npde) should be a column vector.
%
% TauFun:   the function defines the parameter tau. it should have the form
%
%           function tau=TauFun(npde,jmax,t,x,u)
%
%           where (npde,jmax,t,u) are input variable and u contains the values of
%           the physical variables and the mesh at time t. x is a column vector
%           of size jmax and u is a matrix of size npde by jmax: u(i,j), i=1:npde,
%           j=1:jmax.
%
%           default value tau = 0.01 is invoked when TauFun is not defined
%           thru absence or putting [] for this argument.
%
% monitor:  an integer specifying the mesh density function.
%           (also see explanations for variables alpha_type and alpha.)
%           0: for fixed mesh.
%           1: k=0, m=0 (piecewise constant approximation, error measured in L2).
%           2: k=1, m=0 (piecewise linear approximation, error measured in L2).
%           3: k=1, m=1 (piecewise linear approximation, error measured in
%                        H1 semi-norm).
%           4: arclength mesh density function, rho = sqrt(1+(u')^2).
%           5: curvature mesh density function, rho = (1+(u'')^2)^0.25.
%
%           default: monitor=3
%
% reltol, abstol: relative and absolute tolerances used in integration of
%           the ode system.
% MonitorFun: the user defined monitor function. it should have the form
%
%           function rho=MonitorFun(npde,jmax,t,x,u)
%
%           where (npde,jmax,t,x,u) are input variable and u contains the values of
%           the physical variables and the mesh at time t. x is a column vector
%           of size jmax and u is a matrix of size npde by jmax: u(i,j), i=1:npde,
%           j=1:jmax.
%
%           an internal monitor function (specified through (monitor)) is invoked
%           when MonitorFun is not defined thru absence or putting [] for
%           this argument.
%
%   note: when the user provides his own monitor function, the function must return
%         a positive, column vector of size (jmax). in this case, the variable
%         (monitor) should be set to a value other than 0.
% mmpde:    (4, 5, 6, 7) specifies which MMPDE (out of MMPDE4, MMPDE5, MMPDE6,
%           modified MMPDE5) is used. these MMPDEs have the expressions as
%
%           MMPDE5: x_t = 1/(tau) (rho x_xi)_xi,
%           MMPDE4: (rho (x_t)_xi)_xi = -1/(tau) (rho x_xi)_xi,
%           MMPDE6: ((x_t)_xi)_xi = -1/(tau) (rho x_xi)_xi,
%  modified MMPDE5: x_t = 1/(rho*tau) (rho x_xi)_xi,  (default)
%
%           where rho is the monitor function (or more precisely, the mesh density
%           function) and tau, which can be defined as a function of t, mesh,
%           and the physical variables, is a parameter for adjusting the response
%           of the mesh to the changes in rho.
% alpha_type and alpha: these two variables are used to control the definition of
%           the adaptation intensity parameter alpha. alpha_type indicates
%           the type of the definition.
%           alpha_type = 1: integral, solution based definition. in this case,
%                           the value of variable (alpha) is not used.
%           alpha_type = 2: integral, solution based definition with flooring. that is,
%                           alpha := max(alpha, integral definition).
%           alpha_type = 3: alpha := alpha (constant).
%
%           default: alpha_type = 1, alpha = 1
%

%
% Copyright (C) 2010 Weizhang Huang and Robert D. Russell
% all rights reserved.
%
% This program is provided "as is", without warranty of any kind.
% Permission is hereby granted, free of charge, to use this program
% for personal, research, and education purposes. Distribution or use
% of this program for any commercial purpose is permissible
% only by direct arrangement with the copyright owner.
%

% functions contained or called in this file:
%
% 1. res=movfdm_res(t,v,vt,PDE_F,PDE_G,BC_L,BC_R,jmax,npde,TauFun,monitor)
% 2. [SPDY,SPDYP] = jpattern(npde,jmax,mmpde,alpha_type)
% 3. rho = mesh_density_fnct(npde,jmax,x,v,monitor,lsfitting,alpha_type,alpha) 
%                                             (called external function)
% 4. dv = grad_recovery(jmax,x,v)             (called external function)
% 5. rho1 = mesh_density_fnct_smoothing(jmax,rho)      (called external function)
% 6. tau = TauFun_default(npde,jmax,t,x,v)
%

   % set the default values for the parameters

   if nargin < 11
       error('The first 11 arguments must be defined.');
   elseif nargin < 12
      TauFun = @TauFun_default;
      monitor=[];
      reltol=[];
      abstol=[];
      MonitorFun=[];
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 13
      monitor=[];
      reltol=[];
      abstol=[];
      MonitorFun=[];
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 14
      reltol=[];
      abstol=[];
      MonitorFun=[];
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 15
      abstol=[];
      MonitorFun=[];
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 16
      MonitorFun=[];
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 17
	  mmpde=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 18
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 19
	  alpha=[];
   elseif nargin > 19
      error('Two many arguments are specified.');
   end    

   flag_TauFun = 0;
   flag_monitor = 0;
   flag_reltol = 0;
   flag_abstol = 0;
   flag_mmpde = 0;
   flag_alpha_type = 0;
   flag_alpha = 0;
   if isempty(TauFun)
      TauFun = @TauFun_default;
      flag_TauFun = 1;
   end
   if isempty(monitor) | (monitor<0) | (monitor>5)
      monitor=6;
      flag_monitor = 1;
   end
   if isempty(reltol) | (reltol<=0)
      reltol=1e-6;
      flag_reltol = 1;
   end
   if isempty(abstol) | (abstol<=0)
      abstol=1e-4;
      flag_abstol=1;
   end
   if isempty(mmpde) | (mmpde<4) | (mmpde>7)
      mmpde=7;
	  flag_mmpde = 1;
   end
   if isempty(alpha_type) | (alpha_type<1) | (alpha_type>3)
	  alpha_type=1;
	  flag_alpha_type = 1;
   end
   if isempty(alpha) | (alpha<=eps)
      alpha = 1.0;
	  flag_alpha = 1;
   end

% set options

   nn = npde*jmax;
			   
   [SPDY,SPDYP]=jpattern(npde,jmax,mmpde,alpha_type);
   j=round(jmax/2);
   if j < 1, j = 1; end
   j = npde*(j-1)+1;
   opts=odeset('RelTol',reltol,'AbsTol',abstol,...
               'OutputSel',j,'JPattern',{SPDY,SPDYP});
      
% use a helper to find better initial derivatives

   free_v=zeros((npde+1)*jmax,1);
   free_v(nn+1)=1;
   free_v(nn+npde)=1;
   v0 =[reshape(u0,npde*jmax,1)
        x0];
   vt0 =[reshape(vt0,npde*jmax,1)
        x0];
   
   [v0md,vt0]=decic(@movfdm_res,tspan(1),v0,free_v,vt0,[],opts,...
               PDE_F,PDE_G,BC_L,BC_R,...
			   jmax,npde,TauFun,monitor,MonitorFun,mmpde,alpha_type,alpha);
   dd1 = norm(v0md(1:nn)-v0(1:nn),inf);
   dd2 = norm(v0md((nn+1):end)-v0((nn+1):end),inf);
   
% call ode15i

   [t,v]=ode15i(@movfdm_res,tspan,v0md,vt0,opts,...
                PDE_F,PDE_G,BC_L,BC_R,...
				jmax,npde,TauFun,monitor,MonitorFun,mmpde,alpha_type,alpha);
   N=size(t,1);
   x=zeros(jmax,N);
   u=zeros(npde,jmax,N);
   for n=1:N
      u(:,:,n)=reshape(v(n,1:nn),npde,jmax);
      x(:,n)=v(n,(nn+1):end)';
   end
                   
% print messages

if flag_TauFun==1
   fprintf('\n *** The default value tau=0.01 is used.\n');
end
if flag_monitor==1
   fprintf('\n *** The default value monitor=3 is used.\n');
end
if flag_reltol==1
   fprintf('\n *** The default value reltol=1e-6 is used.\n');
end
if flag_abstol==1
   fprintf('\n *** The default value abstol=1e-4 is used.\n');
end
if flag_mmpde==1
   fprintf('\n *** The default mmpde=7 (modified MMPDE5) is used.\n');
end
if flag_alpha_type==1
   fprintf('\n *** The default definition (integral definition) for alpha is used.\n');
end
if flag_alpha==1
   fprintf('\n *** The default alpha = 1.0 is used.\n');
end
if (dd1+dd2)>eps
   fprintf('\n *** Warning: a modified initial solution is used.\n');
     fprintf('     max diff in solution = %e\n', dd1);
     fprintf('         max diff in mesh = %e\n', dd2);
end

% end of movfdm

% --------------------------------------------------------------------------

function res=movfdm_res(t,v,vt,PDE_F,PDE_G,BC_L,BC_R,...
                        jmax,npde,TauFun,monitor,MonitorFun,mmpde,alpha_type,alpha)
%
% this function computes the res function to be called by an ode
% solver. the res function results from the central finite difference
% discretizetion of the (npde) physical pdes and the mesh equation
% and the corresponding boundary conditions.
%
% mmpde = 4, 5, 6, or 7 (modified MMPDE5)
%

   % declare arrays
   
   Res=zeros(npde,jmax);
   Resx=zeros(jmax,1);
   
   u=zeros(npde,1);
   ux=zeros(npde,1);
   ut=zeros(npde,1);
   XX = zeros(1,jmax-1);
   U = zeros(npde,jmax-1);
   UX = zeros(npde,jmax-1);
   UT = zeros(npde,jmax-1);
   G = zeros(npde,jmax-1);   
   rho=zeros(jmax,1);
   nn = npde*jmax;
   
   % extract physical variables and mesh from v and vt

   V = reshape(v(1:nn),npde,jmax);
   Vt = reshape(vt(1:nn),npde,jmax);
   x = v((nn+1):end);
   xt = vt((nn+1):end);

   % compute monitor function and tau
   
   if monitor~=0
      if isempty(MonitorFun) % use internal monitor functions
	     lsfitting = 'no'; % using finite difference approximation
	     rho = mesh_density_fnct(npde,jmax,x,V,monitor,lsfitting,alpha_type,alpha);
      else
         rho = feval(MonitorFun,npde,jmax,t,x,V);
      end      
      % smooth the monitor function   
      for j=1:4
         Resx = rho;
         rho = mesh_density_fnct_smoothing(jmax,Resx);
      end
      % compute tau   
      tau = feval(TauFun,npde,jmax,t,x,V);
      if isempty(tau) | (tau<=0)
         error('tau=TauFun(npde,jmax,t,x,v) should be positive.');
      end
   end

   % for physical pdes: interior points
   
   % for F term
   for j=2:(jmax-1)
	  XX(j) = x(j);
      U(:,j) = V(:,j);
      UX(:,j) = (V(:,j+1)-V(:,j-1))/(x(j+1)-x(j-1));
      UT(:,j) = Vt(:,j)-UX(:,j)*xt(j);
   end
   % Res(:,1) will be reset by boundary conditions
   XX(1) = x(1);
   U(:,1) = V(:,1);
   UX(:,1) = (V(:,2)-V(:,1))/(x(2)-x(1));
   UT(:,1) = Vt(:,1)-UX(:,1)*xt(1);
   Res(:,1:(jmax-1)) = feval(PDE_F,t,XX,U,UX,UT);  
   
   % for G term
   for j=1:(jmax-1)
	  XX(j) = (x(j+1)+x(j))*0.5;
      U(:,j) = (V(:,j+1)+V(:,j))*0.5;
      UX(:,j) = (V(:,j+1)-V(:,j))/(x(j+1)-x(j));
	  UT(:,j) = 0.5*(Vt(:,j+1)+Vt(:,j))-UX(:,j)*0.5*(xt(j+1)+xt(j));
   end
   G = feval(PDE_G,t,XX,U,UX,UT);
   for j=2:(jmax-1)
	  Res(:,j) = Res(:,j) - (G(:,j)-G(:,j-1))*2/(x(j+1)-x(j-1));
   end
   
   % for physical pdes: boundary points

   % at the left end
   j = 1;   
   xx = x(j);
   dx = 1.0/(x(j+1)-x(j));
   dx1 = x(2)-x(1);
   dx2 = x(3)-x(1);
   u = V(:,j);
   % 1st order approximation
%   ux = dx*(V(:,j+1)-V(:,j));
   % 2nd order approximation
   ux = (dx1*dx1*V(:,3)-dx2*dx2*V(:,2)-(dx1*dx1-dx2*dx2)*V(:,1)) ...
        /((dx1-dx2)*dx1*dx2);
   ut = Vt(:,j)-ux*xt(j);
   Res(:,j) = feval(BC_L,t,xx,u,ux,ut);
   
   % at the right end
   j = jmax;   
   xx = x(j);
   dx = 1.0/(x(j)-x(j-1));
   dx1 = x(j)-x(j-1);
   dx2 = x(j)-x(j-2);
   u = V(:,j);
   % 1st order approximation
%   ux = dx*(V(:,jmax)-V(:,jmax-1));
   % 2nd order approximation
   ux = -(dx1*dx1*V(:,j-2)-dx2*dx2*V(:,j-1)-(dx1*dx1-dx2*dx2)*V(:,j)) ...
        /((dx1-dx2)*dx1*dx2);
   ut = Vt(:,j)-ux*xt(j);
   Res(:,j) = feval(BC_R,t,xx,u,ux,ut);
      
   % for mesh equation
   
   if monitor==0 % for fixed meshes

       Resx = xt;
    
   else  % for moving meshes
      
      % for mesh pde: interior points (mmpde5 with balancing factor)
      dxi = jmax-1;
	  if mmpde==4
         for j=2:(jmax-1)
            Resx(j) = (tau*(xt(j+1)-xt(j))+(x(j+1)-x(j)))*(rho(j+1)+rho(j)) ...
              -(tau*(xt(j)-xt(j-1))+(x(j)-x(j-1)))*(rho(j)+rho(j-1));
         end
      elseif mmpde==6
         for j=2:(jmax-1)
            Resx(j) = (tau*(xt(j+1)-xt(j))+(x(j+1)-x(j))*(rho(j+1)+rho(j))) ...
              -(tau*(xt(j)-xt(j-1))+(x(j)-x(j-1))*(rho(j)+rho(j-1)));
         end
      elseif mmpde==5
         for j=2:(jmax-1)
            Resx(j) = tau*xt(j)-0.5*dxi*dxi*((x(j+1)-x(j))*(rho(j+1)+rho(j)) ...
                 -(x(j)-x(j-1))*(rho(j)+rho(j-1)));
         end
      else % for modified MMPDE5 (default)
         for j=2:(jmax-1)
            Resx(j) = tau*xt(j)-0.5*dxi*dxi/rho(j)*((x(j+1)-x(j))*(rho(j+1)+rho(j)) ...
                 -(x(j)-x(j-1))*(rho(j)+rho(j-1)));
         end
      end
              
      % for mesh pde: boundary points (fixed)
   
      Resx(1) = xt(1)-0.0;
      Resx(jmax) = xt(jmax)-0.0;
   end
   
   res = [reshape(Res,nn,1)
          Resx];

% end of movfdm_res

% --------------------------------------------------------------------------

function [SPDY,SPDYP] = jpattern(npde,jmax,mmpde,alpha_type)
%
% this function defines the sparse pattern of the jacobian matrix.
%

   SPDY = sparse((npde+1)*jmax, (npde+1)*jmax);
   SPDYP = sparse((npde+1)*jmax, (npde+1)*jmax);
   nn = npde*jmax;

   for j=1:jmax
      i=(npde*(j-1)+1):(npde*j);
      j1 = max(j-2,1);
      j2 = min(j+2,jmax);	  
      SPDY(i,(npde*(j1-1)+1):(npde*j2)) = ones(npde,(j2-j1+1)*npde);
      SPDY(i,(nn+j1):(nn+j2)) = ones(npde,j2-j1+1);
	  SPDYP(i,i) = ones(npde,npde);
	  SPDYP(i,nn+j) = ones(npde,1);
   end

   if alpha_type==3  % for alpha = constant
      for j=1:jmax
	     j1 = max(j-6,1);
	     j2 = min(j+6,jmax);
	     SPDY(nn+j,(npde*(j1-1)+1):(npde*j2)) = ones(1,(j2-j1+1)*npde);
	     SPDY(nn+j,(nn+j1):(nn+j2)) = ones(1,j2-j1+1);
      end
   else % for alpha based on solution derivatives   
      SPDY((nn+1):end,:) = ones(jmax,nn+jmax);
   end

   e = ones(jmax,1);
   if (mmpde==4) | (mmpde==6) 
      SPDYP((nn+1):end,(nn+1):end) = spdiags([e e e],-1:1,jmax,jmax);
   else
      SPDYP((nn+1):end,(nn+1):end) = spdiags(e,0,jmax,jmax);
   end
   
% end of jpattern

% -----------------------------------------------------

function tau = TauFun_default(npde,jmax,t,x,v)

   tau = 0.01;


   
