function [ ys, check ] = FA_steadystate(ys, exe)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global M_



check = 0;

betta= M_.params(1);
gammaEZ = 2;
psiEZ = 0.5;
thetaEZ = (1-gammaEZ)/(1-1/psiEZ);
nuEZ = 0.3622;

C = log(0.502901560401090);
D = log(0.999999999999998);
F = log(2.68169616988352);
G = log(0.169757100000000);
I = log(0.134617350451695);
In = 0; 
K = log(5.38469424915843);
L = log(0.317031808439407);
N = log(1.34617349036349);
Ne = log(1.33417757620910);
Nn = log(0.0119959141543902);
Pm = log(0.760019198464123);
Q = log(1);
Rk = log(1.01260101027508);
U = log(0.999999997412428);
Y = log(0.807276010852785);
Ym = log(0.807276010852783);
Z = log(3.52845845907945);
delta = log(0.0249999989270949);
eta = log(1.51102093246097);
i = log(1.01010101010101);
infl = log(1.00000000000000);
inflstar = log(1.00000000000000);
nu = log(0.00373977706823913);
phi = log(4.00000021372022);
v =  log(0.611304736762464);
x = log(1.02010101133160);
Lambda = log(1^((1-gammaEZ)/thetaEZ)*1^(1-1/thetaEZ));
X=-Pm;
R = Lambda-log(betta);
a=log(1);
ksi=log(1);
g=log(1);
z=x;

ys = [C D F G I In K L N Ne Nn Pm Q Rk U Y Ym Z delta eta i infl inflstar nu phi v x Lambda X R a ksi g z]';

end

