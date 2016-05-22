function [x_new,nx_new]=static_local_refine(x,u,basex,npdes,nmax,alpha,beta,tolx,bound,imesh,ilim)

% grid_refine computes a new grid, which is equidistributing a monitor function
%
%... input parameters:
%...
%...     z(nz)          : current location of nodes
%...     x(npdes*nz)    : current dependent variable matrix
%...     npdes          : number of PDEs
%...     nzmax          : maximum number of nodes
%...     alpha          : parameter of the monitor function 
%...                     (limit the maximum grid interval, i.e., 
%...                      hmax**2 < tolz**2/alpha)
%...     beta           : parameter of the monitor function
%...                      (limit the excessive clustering of nodes) 
%...     tolz           : reference quantity for the equidistribution 
%...                      of the monitor function
%...     bound >=1      : parameter defining a locally bounded grid
%...     imesh      =0  : monitorfunction f = dqrt( alpha + sum (x_z)^2 )
%...                =1  : monitorfunction f = dsqrt( alpha + max(dabs(x_zz)))
%...     ilim       =0  : no limit on the second spatial derivatives of x
%...                =1  : limit the value of x_zz to beta 
%...
%... output parameters:
%...
%...     nz_new         : new number of nodes
%...     z_new(nz_new)  : new location of nodes
%...


%define size of base grid
nx=length(basex);

% Definition of the monitor function
pp = csapi(x,u);    %use cubic splines to compute spatial derivs
dpp = fnder(pp);
ddpp= fnder(dpp);
ux = fnval(dpp,basex)';  %1st order spatial deriv
uxx = fnval(ddpp,basex)';    %2nd order spatial deriv
       
%depending on imesh, the monitor function is in the form
% mon=sqrt(alpha+||u_x||^2) or mon=||u_xx||
%if imesh = 1, determine max(abs(x_zz))
if imesh == 1
    
    uxxm=abs(uxx); %use the absolute value to determine error
    
    %if ilim = 1, limit the value of max(abs(x_zz))
    if ilim == 1
        index = (uxxm > beta*ones(1,nx));
        uxxm(index) = beta;
    end
end     

%assemble the monitor function
if imesh == 0
    mon=sqrt(alpha+abs(ux).^2);
elseif imesh == 1
    mon=uxxm;
end

% initialise new mesh with the first base grid point
x_new = x(1);

% iterate over all subgrids and define the number of mesh points
for i=1:nx-1
    
    %determine the degree of refinement
    higherErr = max(mon(i),mon(i+1));
    if higherErr <= 1e-1    %check first threshold -> add 0 points
        nsubgrid = 2;
    elseif higherErr <= 1    %second -> add 1 points
        nsubgrid = 3;
    elseif higherErr <= 3    %third -> add 2 points
        nsubgrid = 4;
    elseif higherErr <= 6    %fourth -> add 3 points
        nsubgrid = 5;
    elseif higherErr <= 10    %fifth -> add 4 points
        nsubgrid = 6;
    elseif higherErr <= 15    %sixth -> add 5 points
        nsubgrid = 7;
    elseif higherErr <= 20    %seventh -> add 6 points
        nsubgrid = 8;
    elseif higherErr <= 30    %eighth-> add 7 points
        nsubgrid = 9;
    elseif higherErr <= 40    %nineth -> add 8 points
        nsubgrid = 10;
    elseif higherErr <= 50    %tenth -> add 9 points
        nsubgrid = 11;
    elseif higherErr <= 60    %eleventh -> add 10 points
        nsubgrid = 12;
    else                      %else add 14 points
        nsubgrid = 16;
    end
    
    %create new grid along the interval [x(i), x(i+1)]
    clearvars subgrid
    subgrid(1:nsubgrid)=linspace(basex(i),basex(i+1),nsubgrid);
    
    % update the new mesh
    if length(subgrid) > 1
        x_new=[x_new subgrid(2:end)];
    else
        x_new=[x_new subgrid];
    end
end

nx_new = length(x_new);

