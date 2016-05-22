function rho1 = mesh_density_fnct_smoothing(jmax,rho)
%
% rho1 = adp_fnct_smoothing(jmax,rho)
%
% this function smooths the mesh density function using an averaging scheme.
%
% input variables:
%
% rho:      mesh density function, column vector of size jmax.
%
% output variables:
%
% rho1:     smoothened mesh density function, column vector of size jmax.
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

   rho1 = zeros(jmax,1);
   for j=2:(jmax-1)
      rho1(j) = 0.25*(rho(j-1)+rho(j+1))+0.5*rho(j);
   end
   rho1(1) = 0.5*(rho(1)+rho(2));
   rho1(jmax) = 0.5*(rho(jmax)+rho(jmax-1));

% end of mesh_density_fnct_smoothing
