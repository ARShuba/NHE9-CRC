function out = NHE
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
%-----------------------------------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)

%Ym function components
Mu0=1/(1+kmrgd(1)/u2000)^nu200;
Mu1=(kmrgd(1)/u2000)/(1+kmrgd(1)/u2000)^nu200;
Mu2=(kmrgd(1)/u2000)^2/(1+kmrgd(1)/u2000)^nu200;
Mu3=(kmrgd(1)/u2000)^3/(1+kmrgd(1)/u2000)^nu200;
Mu4=(kmrgd(1)/u2000)^4/(1+kmrgd(1)/u2000)^nu200;
Mu5=(kmrgd(1)/u2000)^5/(1+kmrgd(1)/u2000)^nu200;
Mu6=(kmrgd(1)/u2000)^6/(1+kmrgd(1)/u2000)^nu200;

%Hills functions
Hillszu200=(1+lamdazu200*(kmrgd(3)/z0u200)^nzu200)/(1+(kmrgd(3)/z0u200)^nzu200);
Hillssu200=(1+lamdasu200*(s/s0u200)^nsu200)/(1+(s/s0u200)^nsu200);
Hillszmz=(1+lamdazmz*(kmrgd(3)/z0mz)^nzmz)/(1+(kmrgd(3)/z0mz)^nzmz);
Hillssmz=(1+lamdasmz*(s/s0mz)^nsmz)/(1+(s/s0mz)^nsmz);
Hillsemz = (1+lamdaemz*(kmrgd(4)/e0mz)^nemz)/(1+(kmrgd(4)/e0mz)^nemz);
Hillsze = (1+lamdaze*(kmrgd(3)/z0e)^nze)/(1+(kmrgd(3)/z0e)^nze);
Hillsnz = (1+lamdanz*(kmrgd(5)/n0z)^nnz)/(1+(kmrgd(5)/n0z)^nnz);
Hillsnn = (1+lamdann*(kmrgd(5)/n0n)^nnn)/(1+(kmrgd(5)/n0n)^nnn);
Hillsne = (1+lamdane*(kmrgd(5)/n0e)^nne)/(1+(kmrgd(5)/n0e)^nne);
Hillsnu200=(1+lamdanu200*(kmrgd(5)/n0u200)^nnu200)/(1+(kmrgd(5)/n0u200)^nnu200);
Hillsl7nh = (1+lamdal7nh*(kmrgd(6)/l70nh)^nl7nh)/(1+(kmrgd(6)/l70nh)^nl7nh);
Hillsnhn = (1+lamdanhn*(kmrgd(7)/nh0n)^nnhn)/(1+(kmrgd(7)/nh0n)^nnhn);
Hillsl7l7 = (1+lamdal7l7*(kmrgd(6)/l70l7)^nl7l7)/(1+(kmrgd(6)/l70l7)^nl7l7);
Hillssl7 = (1+lamdasl7*(s/s0l7)^nsl7)/(1+(s/s0l7)^nsl7);

dydt=[gu200*Hillszu200*Hillssu200*Hillsnu200-kmrgd(2)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*kmrgd(1);
gmz*Hillszmz*Hillssmz*Hillsemz*Hillsnz-kmrgd(2)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*kmrgd(2);
gz*kmrgd(2)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*kmrgd(3)
ge*Hillsze*Hillsne-ke*kmrgd(4);
gn*Hillsnn*Hillsnhn-kn*kmrgd(5);
gl7*Hillssl7*Hillsl7l7-kl7*kmrgd(6);
gnh*Hillsl7nh-knh*kmrgd(7);
] ;
    
    % --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(NHE);
y0=[0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
%--------------------------------------------------------------------------
function hess = hessians(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
%--------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,s,ku200,kmz,kz,gu200,gmz,gz,z0u200,z0mz,s0u200,s0mz,u2000,nzu200,nsu200,nzmz,nsmz,nu200,lamdazu200,lamdasu200,lamdazmz,lamdasmz,lamdaemz,e0mz,nemz,lamdaze,z0e,nze,ke,ge,gn,kn,nnz,n0z,lamdanz,nnn,n0n,lamdann,lamdane,n0e,nne,lamdanu200,nnu200,n0u200,lamdal7nh,l70nh,nl7nh,lamdanhn,nh0n,nnhn,kl7,knh,gl7,gnh,lamdal7l7,l70l7,nl7l7,lamdasl7,s0l7,nsl7)
