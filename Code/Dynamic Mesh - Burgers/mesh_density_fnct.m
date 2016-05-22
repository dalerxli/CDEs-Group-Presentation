function rho = mesh_density_fnct(npde,jmax,x,v,monitor,lsfitting,alpha_type,alpha)
%
% function rho = mesh_density_fnct(npde,jmax,t,x,v,monitor,lsfitting,alpha,flooring)
%
% this function computes one of several mesh density functions.
%
% Note that the first 4 arguments must be specified by the user. The other 4
% arguments (monitor,lsfitting,alpha,flooring) are optional.
%
% output variables:
%
% rho:      mesh density function, column vector of size jmax: rho(j), j=1:jmax.
%
% input variables:
%
% npde:     the number of physical PDEs contained in the system to be integrated.
% jmax:     the number of mesh points.
% x:        mesh, column vector of size jmax: x(j), j=1:jmax.
% v:        nodal values of the solution, array of size npde by jmax:
%           v(i,j): i=1:npde, j=1:jmax
% monitor:  an integer specifying the mesh density function.
%           0: for fixed mesh.
%           1: k=0, m=0 (piecewise constant approximation, error measured in L2).
%           2: k=1, m=0 (piecewise linear approximation, error measured in L2).
%           3: k=1, m=1 (piecewise linear approximation, error measured in
%                        H1 semi-norm).
%           4: arclength mesh density function, rho = sqrt(1+(u')^2).
%           5: curvature mesh density function, rho = (1+(u'')^2)^0.25.
%
%           default: monitor=3
% lsfitting: character array or string, indicates if the least squares fitting or
%           finite difference approximation is used to computate approximations
%           to solution deriatives needed in computing the mesh density function.
%           lsfitting = 'yes','y','YES',or 'Y': least squares fitting is used.
%
%           default: lsfitting='yes'
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

%
% functions contained or called in this file:
%
%    dv = grad_recovery(jmax,x,v1)        (called external function)
%

   % set the default values for the parameters
   
   if nargin < 4
       error('The first 4 arguments must be defined.');
   elseif nargin < 5
      monitor=[];
	  lsfitting=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 6
	  lsfitting=[];
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 7
	  alpha_type=[];
	  alpha=[];
   elseif nargin < 8
      alpha=[];
   elseif nargin > 8
      error('Two many arguments are specified.');
   end    

   if isempty(monitor) | (monitor<0) | (monitor>5)
      monitor=3;
   end
   if isempty(lsfitting)
      lsfitting='yes';
   end
   if isempty(alpha_type) | (alpha_type<1) | (alpha_type>3)
	  alpha_type=1;
   end
   if isempty(alpha) | (alpha<=eps)
      alpha = 1.0;
   end
   
   % declare arrays

   rho = zeros(jmax,1);
   v1 = zeros(jmax,1);
   v2 = zeros(jmax,1);
   
   % define k and m

   if monitor==0 % for uniform mesh
      rho=ones(jmax,1);
	  return;
   else
      switch monitor
	     case 1
		    k=0; m=0;
	     case 2
		    k=1; m=0;
	     case 3
		    k=1; m=1;
		 case 4  % arc-length
		    k=0; m=0;
		 case 5 % curvature
		    k=1; m=0;
	     otherwise
		    k=1; m=1;
      end
   end
   
   % compute the (k+1)th order derivatives using either a gradient recovery
   % technique based on the least squares fitting or finite difference
   % approximation. the least suares method is generally fairly robust and
   % accurate but expensive.

   if ~isempty(strmatch(lower(lsfitting),'yes','exact')) ...
      | ~isempty(strmatch(lower(lsfitting),'y','exact'))
      % least squares fitting   
      for l=1:npde
         v2 = v(l,:)';
         for i=1:(k+1)
            v1 = v2;
            v2 = grad_recovery(jmax,x,v1);
         end
         rho = rho + v2.^2;
      end	  
   else % finite difference approximation   
      if k==0 % only 1st derivative is needed.
         for l=1:npde
	        for j=2:(jmax-1)
                v2(j) = (v(l,j+1)-v(l,j-1))/(x(j+1)-x(j-1));
		    end
		    v2(1) = (v(l,2)-v(l,1))/(x(2)-x(1));
		    v2(jmax) = (v(l,jmax)-v(l,jmax-1))/(x(jmax)-x(jmax-1));
            rho = rho + v2.^2;
         end
      else % only 2nd derivative is needed.
         for l=1:npde;
	        for j=2:(jmax-1)
               v2(j)=2.0./(x(j+1)-x(j-1))*((v(l,j+1)-v(l,j))/(x(j+1)-x(j)) ...
                     -(v(l,j)-v(l,j-1))/(x(j)-x(j-1)));
		    end
            v2(1) = 2*((x(2)-x(1))*(v(l,3)-v(l,1))-(x(3)-x(1))*(v(l,2)-v(l,1))) ...
                    /((x(3)-x(1))*(x(2)-x(1))*(x(3)-x(2)));
            v2(jmax) = 2*((x(jmax-1)-x(jmax))*(v(l,jmax-2)-v(l,jmax)) ...
                              -(x(jmax-2)-x(jmax))*(v(l,jmax-1)-v(l,jmax))) ...
                    /((x(jmax-2)-x(jmax))*(x(jmax-1)-x(jmax))*(x(jmax-2)-x(jmax-1)));
            rho = rho + v2.^2;
         end
      end
   end
   
   % compute the arclength mesh density function
   
   if (monitor==4)
      for j=1:jmax
         rho(j)=sqrt(1+rho(j));
	  end
      return;
   end
   
   % compute the curvature mesh density function
   
   if (monitor==5)
      for j=1:jmax
         rho(j)=(1+rho(j))^0.25;
	  end
      return;
   end

   % compute alpha

   gamma = 1.0/(2*(k-m+1)+1);
   Alpha = eps;
   for j=2:jmax
      Alpha = Alpha + 0.5*(rho(j)^gamma+rho(j-1)^gamma)*(x(j)-x(j-1));
   end
   Alpha = (Alpha/(x(jmax)-x(1)))^(1.0/gamma);
   
   if (alpha_type==1) % using the integral definition without flooring
      ;
   elseif alpha_type == 2 % integral definition with flooring
	  if (Alpha<alpha)
	     Alpha=alpha;
      end
   else % alpha = constant
      Alpha=alpha;
   end
   
   % compute rho

   rho = (1+(1/Alpha)*rho).^gamma;
   
%   fprintf('alpha = %e\n', Alpha);
   
% end of mesh_density_fnct

