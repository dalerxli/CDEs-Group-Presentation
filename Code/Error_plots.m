% small script to compare the maximum grid error from the 3 methods used to
% solve burgers' equation

clear
close all

global nz dt

nz = 151;
dt = 0.1;

sprintf('\nMoL results')
[t_MoL,error_MoL] = main_burgers_MoL(nz,dt);
sprintf('\nMoving Mesh results')
[t_moving,error_moving] = main_burgers(0,nz,dt);
sprintf('\nMesh Refinement results')
[t_refine,error_refine] = main_burgers(1,nz,dt);

nz = 151;

sprintf('\nDynamic Method results')
[t_dynam,error_dynam] = sec1_2_burgersFDM(nz,dt);

nz = 151;

clf
Fig5 = figure(5);
set(Fig5,'position',[50 50 900 600])
loglog(t_MoL,error_MoL,'r',t_moving,error_moving,'b',...
    t_refine,error_refine,'m',t_dynam,error_dynam,'k')
xlabel('t');
ylabel('max grid error');
title(sprintf('Adaptation Frequency = %.2f,',dt));
legend('Standard MoL','Moving Mesh','Mesh Refinement','Dynamic Method','location','NorthWest');
grid on

tstep=num2str(dt);
tstep(tstep=='.')=[];
print('-painters','-dpng',strcat(sprintf('images\\burgers_error_nz%d_dt',nz),tstep))
