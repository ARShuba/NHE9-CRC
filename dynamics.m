function dydt = dynamics(t,kmrgd)

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
if t<=200
    s=0;
elseif t>1300
    s=0;
else
    s=410000;
end

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
Hillsl28l7 = (1+lamdal28l7*(kmrgd(8)/l280l7)^nl28l7)/(1+(kmrgd(8)/l280l7)^nl28l7);
Hillsu200l28 = (1+lamdau200l28*(kmrgd(1)/u2000l28)^nu200l28)/(1+(kmrgd(1)/u2000l28)^nu200l28);
Hillsl28l28 = (1+lamdal28l28*(kmrgd(8)/l280l28)^nl28l28)/(1+(kmrgd(8)/l280l28)^nl28l28);
Hillsl7mz = (1+lamdal7mz*(kmrgd(6)/l70mz)^nl7mz)/(1+(kmrgd(6)/l70mz)^nl7mz);
Hillsl7l28 = (1+lamdal7l28*(kmrgd(6)/l70l28)^nl7l28)/(1+(kmrgd(6)/l70l28)^nl7l28);
Hillsl7l7 = (1+lamdal7l7*(kmrgd(6)/l70l7)^nl7l7)/(1+(kmrgd(6)/l70l7)^nl7l7);
Hillssl7 = (1+lamdasl7*(s/s0l7)^nsl7)/(1+(s/s0l7)^nsl7);

dydt=[gu200*Hillszu200*Hillssu200*Hillsnu200-kmrgd(2)*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*kmrgd(1);
gmz*Hillszmz*Hillssmz*Hillsemz*Hillsnz-kmrgd(2)*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*kmrgd(2);
gz*kmrgd(2)*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*kmrgd(3)
ge*Hillsze*Hillsne-ke*kmrgd(4);
gn*Hillsnn*Hillsnhn-kn*kmrgd(5);
gl7*Hillssl7*Hillsl7l7-kl7*kmrgd(6);
gnh*Hillsl7nh-knh*kmrgd(7);
gl28*Hillsu200l28*Hillsl28l28*Hillsl7l28-kl28*kmrgd(8)
] ;