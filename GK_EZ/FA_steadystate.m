function [ ys, check ] = FA_steadystate(ys, exe)
% computes the steady state for the GK_EZ.mod and uses a numerical
% solver (csolve) to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)
global M_
check = 0;
betta= M_.params(1);

gammaEZ = 2;
psiEZ = 0.5;
thetaEZ = (1-gammaEZ)/(1-1/psiEZ);
nuEZ = 0.33;
cte = (0.451610879111627^nuEZ*(1-0.292857916980711)^(1-nuEZ))^(-1);

Y = -0.293404192805828;
Ym = -0.293404192805826;
K = 1.604245986251769;
L = -1.228067712559721;
I = -2.084633064873953;
C = -0.794934356847043;
G = -1.773386687203001;
Q = 0;
Lambda = 0;
Rk = 0.012522278257187;
R = 0.010050335853501;
N = 0.217951571701830;
Ne = 0.209000509732286;
Nn = -4.502503718808295; 
nu = -5.588729276755017; 
eta = 0.412785536577148;
phi = 1.386294414549939;
z = 0.019901653110221; 
x = 0.019901653110221; 
Pm = -0.274411584883814;
U  = -2.587572048868168*10^(-09);
X  = 0.274411584883814;
D  = 0;
F  = 0.907134949035708;
Z  = 1.181546533919522;
i  = 0.010050335853501;
delta = -3.688879497030137; 
In = 0; 
a  = 0;
ksi = 0;
g  = 0;
infl  = 0;
inflstar = 0;  
v  = 1;
u  = 1;
ev = 1;

ys = [Y Ym K L I C G Q Lambda Rk R N Ne Nn nu eta phi z x Pm U X D F Z i ...
    delta In a ksi g infl inflstar v u ev]';


