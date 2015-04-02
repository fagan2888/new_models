function [ ys, check ] = FA_steadystate(ys, exe)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global M_

check = 0;

betta = M_.params(1);
nuEZ = M_.params(34);

initial= getInitial();
[SS_csolve, rc]=csolve(@FA_stst_csolve, initial,[],1e-12,10000);

if (rc==0)
C = SS_csolve(1);
D = SS_csolve(2);
F = SS_csolve(3);
G = SS_csolve(4);
I = SS_csolve(5);
In = SS_csolve(6); 
K = SS_csolve(7);
L = SS_csolve(8);
N = SS_csolve(9);
Ne = SS_csolve(10);
Nn = SS_csolve(11);
Pm = SS_csolve(12);
Q = SS_csolve(13);
Rk = SS_csolve(14);
U = SS_csolve(15);
Y = SS_csolve(16);
Ym = SS_csolve(17);
Z = SS_csolve(18);
delta = SS_csolve(19);
eta = SS_csolve(20);
i = SS_csolve(21);
infl = SS_csolve(22);
inflstar = SS_csolve(23);
nu = SS_csolve(24);
phi = SS_csolve(25);
x = SS_csolve(25);
varrho = nuEZ * C^(1-nuEZ) * (1-L)^(1-nuEZ);

Lambda = 1;
X=-Pm;
R = Lambda-log(betta);
a=log(1);
ksi=log(1);
g=log(1);
z=x;
elseif rc == 1 | rc == 3
 msg = 'Error occurred steady state not reached. no solution despite extremely fine adjustments in step length (very likely a numerical problem, or a discontinuity';
    error(msg);
elseif rc == 4
    msg = 'Error occurred steady state not found, maximum of Iterations reached';
    error(msg);
end;

ys = [C D F G I In K L N Ne Nn Pm Q varrho Rk U Y Ym Z delta eta i infl inflstar nu phi x Lambda X R a ksi g z]';

ys = log(ys);
ys(6) = exp(ys(6));
end

