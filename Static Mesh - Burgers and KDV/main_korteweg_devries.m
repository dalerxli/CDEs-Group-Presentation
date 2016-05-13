close all;
clear;
% Start a stopwatch timer
tic

%use local_refine function or agereg function?
use_local_refine=0;

% Set global variables
global s
global z0 zL nz D1
global nsteps maxsteps tprint tflag tout
% Spatial grid
z0 = -30.0;
zL = 70.0;
nz = 51;
nzout = nz;
dz = (zL-z0)/(nz-1);
z = [z0:dz:zL]';
basez=z;
error=[];

% Initial conditions
s = 0.5;
x = kdv3_exact(z, 0);
% parameters of the adaptive grid
npdes = 1;
nzmax = 1001;
alpha = 0;
beta = 100;
tolz = 0.005;
bound= 1.1;
%imesh= 0;
ilim = 0;

% refine the initial grid
if use_local_refine == 0
    imesh= 0;
    [z_new, nz_new, ier, tolz_new] = agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
    tolz = tolz_new;
else
    imesh= 1;
    [z_new, nz_new] = static_local_refine(z,x,basez,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
end

% interpolate the dependent variables
x_new = spline(z,x,z_new);
x = x_new';
z = z_new';
nz = nz_new;

% differentiation matrix
D1 = three_point_centered_D1(z);

% call to ODE solver
t0 = 0;
tf = 100;
dt = 5;
yout = x;
zout = z;
nzout= [nzout; nz];
tout = t0;

% solver to stop after this many steps:
maxsteps = 10;
% initial situation
figure(1)
subplot('position',[0.1, 0.3, 0.8, 0.6]);
plot(z,x,'.r-');
ylabel('u(x,t)');
axis([-30, 70, -0.02, 0.27]);
grid on
hold off
subplot('position',[0.1, 0.08, 0.8, 0.17]);
plot(z,t0*ones(nz,1),'.b');
ylabel( 't' );
xlabel( 'x' );
axis([-30, 70, 0, tf]);
grid on
hold on
%
tk = t0;
tspan= [t0, tf];
tprint = dt;
while tk <= tf - 1.e-5
% initialize step counter
	nsteps= 0;

	% do the integration for maxsteps steps in a loop % until t becomes larger than tf
	options = odeset('RelTol',1e-3,'AbsTol',1e-6);
	options = odeset(options,'Events',@count_steps);
	options = odeset(options,'JPattern',jpattern(nz));
	[t,y,te,ye,ie] = ode23s(@kdv3_adaptive_pde,tspan, ...
	 	x,options);
	%
	tk = t(end);
	tspan = [tk, tf];
	x = [];
	x = y(end,:);
    
	% refine the grid
	if use_local_refine == 0
        [z_new, nz_new, ier, tolz_new] = agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
        tolz = tolz_new;
    else
        [z_new, nz_new] = static_local_refine(z,x,basez,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
    end
    
    % interpolate the dependent variables
	x_new = spline(z,x,z_new);
	x =	x_new';
	yout = [yout; x];
    
	z = z_new';
	nz = nz_new;
	zout = [zout; z];
	nzout= [nzout; nz];
	tout = [tout ; tk];

	% plot intermediate results
	if tflag >= 0
		figure(1)
		subplot('position', [0.1, 0.3, 0.8, 0.6]);
        pause(.1)
        plot(z,x,'.r-','markersize', 10);
        axis([-30, 70, -0.02, 0.27]);
        grid on
        hold on
        
	 	clearvars uexact
		uexact = kdv3_exact(z,tk);
		plot(z,uexact(1:length(z)),'b')
        grid on
        hold off
		subplot('position', [0.1, 0.08, 0.8, 0.17])
		plot(z,tk*ones(nz,1),'.b')
		tprint = tprint + dt;
    end
    
    error=[error max(abs(uexact-x))];
	% compute a new differentiation matrix
	D1 = three_point_centered_D1(z);
end
%output highest error
fprintf('worst error of all mesh points and all time steps = %6.4f\n',max(error));
% read the stopwatch timer
tcpu = toc;
nav = sum(nzout)/length(nzout);
fprintf('time = %6.4f, \naverage # mesh points = %8.4f \n\n',tcpu,nav);
