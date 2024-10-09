function [x,v,s,h,f] = NHE_Bifur

curdir = pwd;
init;
cd(curdir);

opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',200000);
opt=contset(opt,'MinStepsize',1);
opt=contset(opt,'MaxStepsize',100);
opt=contset(opt,'Eigenvalues',1);

% Degradation rate:
ku34 = 0.05;   kms = 0.5;   ks = 0.125;  ku200 = 0.05;   kmz = 0.5;   kz = 0.1;  ke=0.1;
kn=0.1; kl7=0.05; knh=0.1; kl28=0.1;
% Transcription rate:
gu34 = 1350;   gms = 90;   gs = 100; gu200 = 2100;   gmz = 11;   gz = 100; ge=5000; gn=80000;
gl7=200; gnh=80000; gl28=200;
% Hills function threshold :
s0u34 = 300000;   s0ms = 200000;   z0u34 = 600000; u034 = 10000;  e0mz=20000; z0e=100000;
I0ms=50000; z0u200 = 220000;   z0mz = 25000;   s0u200 = 180000;   s0mz = 180000; u2000 = 10000;
n0z=800000; n0n=100000; n0e=200000; n0u200=500000; l70nh=200000; nh0n=180000; u2000l28=25000;
l280l28=300000; l70l28=25000; l70l7=1200; l280l7=500000; l70mz=50000; s0l7=180000;
% Cooperativity:
nsu34 = 1;    nsms = 1;   nsmz = 2;   nu34 = 2; nI = 2; nze=2; nnu200=4;
nzu200 = 3;   nsu200 = 2;   nzmz = 2;   nu200 = 6; nzu34=2; nemz=2; nnz=2; nnn=3; nne=4;
nnhn=2; nl7nh=2; nu200l28=2; nl28l28=7; nl7l28=1; nl7l7=3; nl28l7=2; nl7mz=2; nsl7=2;
% fold change
lamdasu34 =0.1;   lamdazu34 = 0.2;  lamdasms = 0.1;   lamdaIms = 10; lamdaze =0.1;
lamdazu200 =0.1;   lamdasu200 = 0.1;  lamdazmz = 7.5;   lamdasmz = 10; lamdaemz=0.8;
lamdanz=2; lamdann=7; lamdane=8; lamdanu200=4; lamdal7nh=0.6;lamdanhn=2;
lamdau200l28=0.5; lamdal28l28=3; lamdal7l28=0.1; lamdal7l7=11; lamdal28l7=0.1; lamdal7mz=0.3;
lamdasl7=0.6;
% external signal
s = 750000;

ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@NHE);
tspan = 0:100:50000;

% initial condition
x_start = [33554.833280 56.500562 0 0 0 0 0];

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7),tspan,x_start);
x_init = x_time(end,:)';

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@NHE,x_init,[s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);