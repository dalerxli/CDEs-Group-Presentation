function [t,error] = main_korteweg_devries()

close all;
clear;
% Start a stopwatch timer
tic

%use local_refine function(=1) or agereg function(=0)?
use_local_refine=0;

%set default axes properties
set(0,'DefaultAxesFontSize', 14,'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',1.5);
%set default legend properties
set(0,'DefaultLegendFontSize',12,'DefaultLegendLocation','NorthWest','DefaultLegendFontSizeMode','manual');
%set default line properties
set(0,'DefaultLineLineWidth',1.5);
%set default grid properties
set(0,'DefaultAxesGridLineStyle','-','DefaultAxesMinorGridLineStyle',':');

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
dt = 2;
yout = x;
zout = z;
nzout= [nzout; nz];
tout = t0;

% solver to stop after this many steps:
maxsteps = 10;
% initial situation
Fig1=figure(1);
set(Fig1, 'Position', [100 100 1000 550])
subplot('position',[0.1, 0.3, 0.85, 0.65]);
plot(z,x,'.b-','MarkerSize',15);
ylabel('u(x,t)');
axis([-30, 70, -0.02, 0.27]);
grid on
hold on
subplot('position',[0.1, 0.1, 0.85, 0.17]);
plot(z,t0*ones(nz,1),'.b');
ylabel( 't' );
xlabel( 'x' );
axis([-30, 70, 0, tf]);
grid on
hold on

uexact = kdv3_exact(z,t0); %exact solution
subplot('position', [0.1, 0.3, 0.85, 0.65]);
plot(z,uexact,'r','LineWidth',2);
axis([-30, 70, -0.02, 0.27]);
legend('Approximate Solution','Exact Solution');
set(gca,'XTickLabel','');
hold off

%save plot
% plot_count=0;
% if use_local_refine == 0
%     print('-painters','-dpng',sprintf('imagesKDV\\static_moving_KDV_dt%d_%d',dt,plot_count))
% else
%     print('-painters','-dpng',sprintf('imagesKDV\\static_adapt_KDV_dt%d_%d',dt,plot_count))
% end

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
		subplot('position', [0.1, 0.3, 0.85, 0.65]);
        pause(.1)
        plot(z,x,'.b-','MarkerSize',15);
        grid on
        hold on
        
	 	clearvars uexact
		uexact = kdv3_exact(z,tk);
		plot(z,uexact(1:length(z)),'r','LineWidth',2)
        ylabel('u(x,t)');
        axis([-30, 70, -0.02, 0.27]);
        legend('Approximate Solution','Exact Solution');
        set(gca,'XTickLabel','');
        grid on
        hold off
        
		subplot('position', [0.1, 0.1, 0.85, 0.17])
		plot(z,tk*ones(nz,1),'.b')
        
        tprint = tprint + dt;
        
        %save plot
%         plot_count = plot_count + 1;
%         if use_local_refine == 0
%             print('-painters','-dpng',sprintf('imagesKDV\\static_moving_KDV_dt%d_%d',dt,plot_count))
%         else
%             print('-painters','-dpng',sprintf('imagesKDV\\static_adapt_KDV_dt%d_%d',dt,plot_count))
%         end
% 		
    end
    error=[error max(abs(uexact-x))];
	% compute a new differentiation matrix
	D1 = three_point_centered_D1(z);
end
% plot and output highest error
figure(2)
loglog(tout(2:end),error,'.b-','markersize',15)
xlabel('t');
ylabel('max grid error');
grid on
% save error plot
if use_local_refine == 0
    print('-painters','-dpng',sprintf('imagesKDV\\static_moving_KDV_error_dt%d',dt))
else
    print('-painters','-dpng',sprintf('imagesKDV\\static_adapt_KDV_error_dt%d',dt))
end

fprintf('maximum grid error across all time steps = %6.4f\n',max(error));
% read the stopwatch timer
tcpu = toc;
nav = sum(nzout)/length(nzout);
fprintf('time = %6.4f, \naverage # mesh points = %8.4f \n\n',tcpu,nav);


t = tout(2:end);
end
