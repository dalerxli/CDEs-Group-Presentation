close all;
clear;
% Start a stopwatch timer
tic

%use local_refine function or agereg function?
use_local_refine=0;

% Set global variables
global mu eps
global z0 zL nz D1
global nsteps maxsteps tprint tflag tout
% Spatial grid
z0 = 0;
zL = 1;
nz = 51; %number of mesh points
nzout = nz; %vector to store all grid sizes used
dz = (zL-z0)/(nz-1); %initial spatital grid step
z = [z0:dz:zL]'; %spatial grid
basez=z;
error=[];

% Initial conditions
mu = 0.001;
for i=1:nz
    x(i) = burgers_exact(z(i),0);
end;

% parameters of the adaptive grid (for agereg function)
npdes = 1;  %number of pdes
nzmax = 1001;   %max grid size
alpha = 0;      %monitor function parameter
beta = 100;     %monitor function parameter
tolz = 0.005;   %equidistribution parameter
bound= 1.1;     %defines a locally bounded grid
imesh= 1;   %(binary) specify monitor funciton
ilim = 1;   %(binary) whether to limit u_xx

% refine the initial grid
if use_local_refine == 0
    [z_new, nz_new, ier, tolz_new] = agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
    tolz = tolz_new;
else
    [z_new, nz_new] = static_local_refine(z,x,basez,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
end

% interpolate the dependent variables
x_new = spline(z,x,z_new);
x = x_new';
z = z_new';
nz = nz_new;

% differentiation matrix
D1 = three_point_centered_D1(z); %use 3 point central FDM

% call to ODE solver
t0 = 0;     %initial time
tf = 1;     %end time
dt = 0.1;     %time step
yout = x;
zout = z;
nzout= [nzout; nz];
tout = t0;

% solver to stop after this many steps:
maxsteps = 5;
% initial situation
Fig1=figure(1);
set(Fig1,'Units','Normalized','OuterPosition',[0 0 1 1]);
subplot('position',[0.1, 0.3, 0.8, 0.6]);
plot(z,x,'.r-','markersize', 10);
ylabel('u(x,t)');
axis([0, 1, 0, 1.1]);
grid on
hold on
subplot('position',[0.1, 0.08, 0.8, 0.17]);
plot(z,t0*ones(nz,1),'.b');
ylabel( 't' );
xlabel( 'x' );
axis([0, 1, 0, tf]);
grid on
hold on

for i=1:nz
    uexact(i) = burgers_exact(z(i),t0); %exact solution
end
subplot('position', [0.1, 0.3, 0.8, 0.6]);
plot(z,uexact,'b');
axis([0, 1, 0, 1.1]);
hold off

%save plot
plot_count=0;
if use_local_refine == 0
    tstep=num2str(dt);
    tstep(tstep=='.')=[];
    print('-painters','-dpng',sprintf(strcat('imagesburg\\static_moving_burgers_dt',...
        tstep,'_%d'),plot_count))
else
    tstep=num2str(dt);
    tstep(tstep=='.')=[];
    print('-painters','-dpng',sprintf(strcat('imagesburg\\static_adapt_burgers_dt',...
        tstep,'_%d'),plot_count))
end


tk = t0;
tspan= [t0, tf];
tprint = dt;
eps = 0.001;

while tk <= tf - 1.e-5
% initialize step counter
	nsteps= 0;

	% do the integration for maxsteps steps in a loop % until t becomes larger than tf
	options = odeset('RelTol',1e-5,'AbsTol',1e-5);
	options = odeset(options,'Events',@count_steps);
	options = odeset(options,'JPattern',jpattern(nz));
	[t,y,te,ye,ie] = ode15s(@burgers_adaptive_pde,tspan, ...
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
	%x_new = spline(z,[1 x(2:end-1) 0.1],z_new);
    x_new = spline(z,x,z_new);
	x =	x_new';
	yout = [yout; x];
    
	z = z_new';
	nz = nz_new;
	%zout = [zout; z];
	nzout= [nzout; nz];
	tout = [tout ; tk];

	% plot intermediate results
	if tflag >= 0
		figure(1)
		subplot('position', [0.1, 0.3, 0.8, 0.6]);
        pause(.1)
	 	plot(z,x,'.r-','markersize', 10);
        axis([0, 1, 0, 1.1]);
        grid on
        hold on

        clearvars uexact
        for i=1:nz
            uexact(i) = burgers_exact(z(i),tk); %exact solution
        end;
		plot(z,uexact(1:length(z)),'b')
        ylabel('u(x,t)');
        grid on
        hold off
        axis([0, 1, 0, 1.1]);
		subplot('position', [0.1, 0.08, 0.8, 0.17])
		plot(z,tk*ones(nz,1),'.b')
		tprint = tprint + dt;
        
        %save plot
        plot_count = plot_count + 1;
        if use_local_refine == 0
            tstep=num2str(dt);
            tstep(tstep=='.')=[];
            print('-painters','-dpng',sprintf(strcat('imagesburg\\static_moving_burgers_dt',...
                tstep,'_%d'),plot_count))
        else
            tstep=num2str(dt);
            tstep(tstep=='.')=[];
            print('-painters','-dpng',sprintf(strcat('imagesburg\\static_adapt_burgers_dt',...
                tstep,'_%d'),plot_count))
        end
		
        
    end
    
    error=[error max(abs(uexact-x'))];
	% compute a new differentiation matrix
	D1 = three_point_centered_D1(z);
end
%output highest error
fprintf('worst error of all mesh points and all time steps = %6.4f\n',max(error));
% read the stopwatch timer
tcpu = toc; %output computation ttime
nav = sum(nzout)/length(nzout); %output average number of mesh points
fprintf('time = %6.4f, \naverage # mesh points = %8.4f \n\n',tcpu,nav);
