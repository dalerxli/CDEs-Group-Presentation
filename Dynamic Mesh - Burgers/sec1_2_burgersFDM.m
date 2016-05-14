function sec1_2_burgersFDM
%
% example driver for Burgers' equation in section 1.2.
% it calls movfdm().
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


clear
clf
cpu0=clock;

   % job = 1 for solution
   % job = 2 for mesh trajectories

   job = 1;
    
   %define number of mesh points
   jmax = 51;
   npde = 1;
   nn = npde*jmax;

   x=zeros(jmax,1);
   u=zeros(npde,jmax);

% define initial solution
   
   x=linspace(0,1,jmax)';
   u(1,:)=(sin(2*pi*x)+0.5*sin(pi*x))';
   xt=zeros(jmax,1);
   ut=zeros(npde,jmax);
   
% Initial conditions
   global epsilon
   epsilon = 0.001;
   for i=1:jmax
      u(1,i) = burgers_exact(x(i,1),0);
   end;
   

   
% define time parameters
t0 = 0;
tf = 1;
dt = 0.1;
% initial situation
Fig1=figure(1);
set(Fig1,'Units','Normalized','OuterPosition',[0 0 1 1]);

subplot('position',[0.1, 0.3, 0.8, 0.6]);
plot(x(:,1),u(1,:),'.r-','markersize', 10);
ylabel('u(x,t)');
axis([0, 1, 0, 1.1]);
grid on
hold on
subplot('position',[0.1, 0.08, 0.8, 0.17]);
plot(x(:,1),t0*ones(jmax,1),'.b');
ylabel( 't' );
xlabel( 'x' );
axis([0, 1, 0, tf]);
grid on
hold on

for i=1:jmax
    uexact(i) = burgers_exact(x(i,1),t0);
end;

subplot('position', [0.1, 0.3, 0.8, 0.6]);
plot(x(:,1),uexact,'b');
axis([0, 1, 0, 1.1]);
hold off

%save plot
print('-painters','-dpng',sprintf('images\\dynam_burgers_%d',0))
%saveas(gcf,sprintf('dynam_burgers%d',0),'pdf');

   
% call the moving mesh function
% monitor = 0 (fixed mesh), 1, 2, 3, 4, or 5
% mmpde = 4, 5, 6, 7 (mmpde = 7 ==> modified MMPDE5)
% alpha_type = 1, 2, or 3: (integral def, integral def with flooring, or alpha = constant)

   monitor = 3;
   mmpde = 7;
   alpha_type = 2;
   alpha = 1.0;
   reltol=1e-4;
   abstol=1e-5;
    
   error=[];
   
   if job==1 % for solution

      tspan = t0:dt:tf;
      
      %solve PDE
      [t,X,U]=movfdm(npde,jmax,tspan,x,xt,u,ut,@PDE_F,@PDE_G,@BC_L,@BC_R,...
                   [],monitor,reltol,abstol,[],mmpde,alpha_type,alpha);
               
      N=size(t,1);

      for n=1:N
         figure(1)
         subplot('position', [0.1, 0.3, 0.8, 0.6]);
         pause(.2)
         
         u=U(:,:,n);
         x=X(:,n);
         

         plot(x,u(1,:)','.r-','markersize', 10);%,'LineWidth',2);
         ylabel('u(x,t)');
         axis([0, 1, 0, 1.1]);
         grid on
         hold on
         
         for i=1:jmax
            uexact(i) = burgers_exact(x(i,1),t(n)); %exact solution
         end;
		 plot(x(:,1),uexact,'b')
         grid on
         hold off
         axis([0, 1, 0, 1.1]);
		 subplot('position', [0.1, 0.08, 0.8, 0.17])
		 plot(x(:,1),t(n)*ones(jmax,1),'.b')
         
         %save plot
         print('-painters','-dpng',sprintf('images\\dynam_burgers_%d',n))
         %saveas(gcf,sprintf('dynam_burgers%d',n),'pdf');
         
         %find max grid error
         error(n)=max(abs(uexact-u(1,:)));
      end

   else % for trajectories
   
      tspan = [0 1.0];
      [t,X,U]=movfdm(npde,jmax,tspan,x,xt,u,ut,@PDE_F,@PDE_G,@BC_L,@BC_R,...
                   [],monitor,reltol,abstol,[],mmpde,alpha_type,alpha);
      plot(X',t);
	  axis([0 1 0 1]);
      xlabel('x')
      ylabel('t')
   end
%output hightest error
fprintf('worst error of all mesh points and all time steps = %6.4f\n',max(error));
%output cputime
fprintf('cpu time used = %e \n\n', etime(clock,cpu0));

% -----------------------------------------------------

function f=PDE_F(t,x,u,ux,ut)
  f(1,:) = ut(1,:);
end
  
function f=PDE_G(t,x,u,ux,ut)
  f(1,:) = epsilon*ux(1,:)-u(1,:).*u(1,:)*0.5;
end
  
function f=BC_L(t,x,u,ux,ut)
  f(1) = ut(1)-0;
end
  
function f=BC_R(t,x,u,ux,ut)
  f(1) = ut(1)-0;
end
  
  
end