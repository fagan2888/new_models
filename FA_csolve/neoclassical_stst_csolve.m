% function yf=neoclassical_stst_csolve(x)
% 
% BETTA = 0.99; %discount rate
% DELTA = .025; %depreciation rate
% ALFA  = 0.33; %capital share
% RHO   = 0.95; %persistence of technology shock
% SIG   = 2;    %intertemporal elasticity of substitution
function yf=neoclassical_stst(x,BETTA,SIG,ALFA,DELTA)

[rows,cols]=size(x);
j=1;
while j<=cols;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% variables in the system:
c    = x(1,j); cp    = x(1,j); 
k    = x(2,j); kp    = x(2,j);
inve = x(3,j); invep = x(3,j);
y    = x(4,j); yp    = x(4,j);
a = 1; ap = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM:
yf(1,j) =   c^(-SIG) - BETTA * cp^(-SIG) * (ap * ALFA * kp^(ALFA-1) + 1 - DELTA);
yf(2,j) =   c + kp - (1-DELTA) * k - a * k^ALFA;
yf(3,j) =   y - a * k^ALFA;
yf(4,j) =   inve - (kp - (1-DELTA)*k);
% y(5,j) =   log(ap) - RHO * log(a);

j=j+1;
end;

% call from Matlab Command Window (or from within some m-file by defining a starting vector x0
% 
% x0 = ones(4,1)*0.9;
%
% and calling the fsolve-function on the system defined by yf:
%
% [SS_VALUES,rc]=csolve(@neoclassical_stst_csolve,x0,[],1e-12,1000), or:
% [SS_VALUES,rc]=csolve(@neoclassical_stst_csolve,x0,[],1e-12,1000,BETTA,SIG,ALFA,DELTA) 
% if you pass through the parameters to neoclassical_stst_fsolve.m
%
% type help csolve to find out how to modify the options
