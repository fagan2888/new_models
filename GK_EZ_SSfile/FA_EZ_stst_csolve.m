% Author: Max Lampe
% Last Change: 17.03.2015
% to Do: 
% - function for csolve steady state
% - trnsformed variables back into levels
% - reduce system by extra variables: Keff, w, VMPK, prem, Welf  
% - Loading parameters from M_.params makes it compatible with dynare and
% automatic calling

function yf=FA_EZ_stst_csolve(n)
[rows,cols] = size(n);
j=1;

global M_

% parameters in the system:

betta = M_.params(1);
sig =  M_.params(2);
varphi =  M_.params(3);
zetta =  M_.params(4);
theta =  M_.params(5);
alfa =  M_.params(6);
G_over_Y =  M_.params(7);
eta_i =  M_.params(8);
epsilon =  M_.params(9);
gam =  M_.params(10);
gam_P =  M_.params(11);
kappa_pi =  M_.params(12);
kappa_y = M_.params(13);
rho_i  =  M_.params(14);
rho_ksi =  M_.params(15);
sigma_ksi =  M_.params(16);
rho_a =  M_.params(17);
sigma_a =  M_.params(18);
rho_g =  M_.params(19);
sigma_g =  M_.params(20);
sigma_Ne =  M_.params(21);
sigma_i = M_.params(22);
rho_shock_psi = M_.params(23);
sigma_psi = M_.params(24);
kappa = M_.params(25);
tau = M_.params(26);
omega = M_.params(27);
lambda = M_.params(28);
chi = M_.params(29);
b = M_.params(30);
delta_c = M_.params(31);
G_ss       =   M_.params(32);
I_ss       =   M_.params(33);

nuEZ = M_.params(34);
psiEZ = M_.params(35);
gammaEZ = M_.params(36); 
thetaEZ = M_.params(37);

while j<=cols;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


%%%%%%%%%%%%% reducing the system of equation %%%%%%%%%%%%% 


% variables in the system:
count=1;
C = n(count,j); count = count+1;
D = n(count,j); count = count+1;
F = n(count,j); count = count+1;
G = n(count,j); count = count+1;
I = n(count,j); count = count+1;
In = n(count,j); count = count+1;
K = n(count,j); count = count+1;
L = n(count,j); count = count+1;
N = n(count,j); count = count+1;
Ne = n(count,j); count = count+1;
Nn = n(count,j); count = count+1;
Pm = n(count,j); count = count+1;
Q = n(count,j); count = count+1;
Rk = n(count,j); count = count+1;
U = n(count,j); count = count+1;
Y = n(count,j); count = count+1;
Ym = n(count,j); count = count+1;
Z = n(count,j); count = count+1;
delta = n(count,j); count = count+1;
eta = n(count,j); count = count+1;
i = n(count,j); count = count+1;
infl = n(count,j); count = count+1;
inflstar = n(count,j); count = count+1;
nu = n(count,j); count = count+1;
phi = n(count,j); count = count+1;
v = n(count,j); count = count+1;
x = n(count,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=1;
% Stochastic discount rate
Lambda = 1^((1-gammaEZ)/thetaEZ)*1^(1-1/thetaEZ);

% Markup
X=1/Pm;
% Euler Eequation for Capital
R = 1 /(Lambda*betta);
a=1;
ksi=1;
g=1;
z=x;

%// 0.1 Utility
u = C^nuEZ*(1-L)^(1-nuEZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM:

% EZ Zin Preferences
%1% // 0.3 Value Function
yf(count,j) = ((1-betta)*u^((1-gammaEZ)/thetaEZ) + betta*((v^(1-gammaEZ))^(1/thetaEZ)))^(thetaEZ/(1-gammaEZ))- v;
count = count +1;

%2% // 0.4 Static Leisure-Consumption
yf(count,j) = Pm*(1-alfa)*Y/L - (1-nuEZ)/nuEZ*C/(1-L);
count = count +1;

% //Financial Intermediaries
%5%  //5. Value of banks' capital
 yf(count,j) = (1-theta)*betta*Lambda*(Rk-R)+betta*Lambda*theta*x*nu-nu;
count = count +1;

% //6. Value of banks' net wealth
yf(count,j) = (1-theta) + betta * Lambda * theta * z * eta - eta;
count = count +1;

% //7. Optimal leverage
yf(count,j) = eta/(lambda-nu) - phi;
count = count +1;

%6% //8. Growth rate of banks' capital
yf(count,j) = (Rk-R) * phi + R - z;
count = count +1;

% //9. Growth rate of banks' net wealth
% yf(count,j) = z-x;
% count = count +1;

% //Aggregate capital, net worth
%7% //10. Aggregate capital
yf(count,j) = phi * N - Q * K;
count = count +1;

%8% //11. Banks' net worth
yf(count,j) = Ne+Nn-N;
count = count +1;

%9% //12. Existing banks' net worth accumulation
yf(count,j) = theta*z*N-Ne;
count = count +1;

% //13. New banks' net worth
yf(count,j) = omega*Q*ksi*K-Nn;
count = count +1;

% //Final goods producer
%10% //14. Return to capital
yf(count,j) = (Pm * alfa * Ym/K + ksi * (Q-delta) )/ Q - Rk;
count = count +1;

%11% //15. Production function
yf(count,j) = a*(ksi*U*K)^alfa*L^(1-alfa) - Ym;
count = count +1;

% //Capital Goods Producer
%12% //16. Optimal investment decision
yf(count,j) = 1+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2 ...
+eta_i*((In+I_ss)/(In+I_ss)-1)*(In+I_ss)/(In+I_ss)-betta*Lambda*eta_i...
*((In+I_ss)/(In+I_ss)-1)*((In+I_ss)/(In+I_ss))^2 - Q;
count = count +1;

%13% //17. Depreciation rate
yf(count,j) = delta_c+b/(1+zetta)*U^(1+zetta) - delta;
count = count +1;

%14% //18. Optimal capacity utilization rate
yf(count,j) = b*U^zetta*ksi*K - Pm*alfa*Ym/U;
count = count +1;

%15% //19. Net investment
yf(count,j) = I-delta*ksi*K - In;
count = count +1;

%16% //20. Capital accumulation equation
yf(count,j) = ksi*K+In-K; 
count = count +1;

%17% //21. Government consumption
yf(count,j) = G_ss*g-G;
count = count +1;

% //Equilibrium
%18% //22. Aggregate resource constraint
yf(count,j) = C+G+I+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2*(In+I_ss) - Y;
count = count +1;

%19% //23. Wholesale, retail output
yf(count,j) = Y*D-Ym;
count = count +1;

%20% //24. Price dispersion
yf(count,j) = gam*D*infl^(-gam_P*epsilon)*infl^epsilon ...
+(1-gam)*((1-gam*infl^(gam_P*(1-gam))*infl^(gam-1))/(1-gam))^(-epsilon/(1-gam)) - D;
count = count +1;

% //25. Markup
% yf(count,j) = 1/Pm-X;
% count = count +1;

%21% //26. Optimal price choice
yf(count,j) = Y*Pm+betta*gam*Lambda*infl^epsilon*(infl)^(-epsilon*gam_P)*F- F;
count = count +1;

%22% //27.
yf(count,j) = Y+betta*gam*Lambda*infl^(epsilon-1)*infl^(gam_P*(1-epsilon))*Z-Z;
count = count +1;

%23% //28. Optimal price choice
yf(count,j) = epsilon/(epsilon-1)*F/Z*infl-inflstar;
count = count +1;

%24% //29. Price index
yf(count,j) = gam*infl^(gam_P*(1-epsilon))+(1-gam)*(inflstar)^(1-epsilon) - infl^(1-epsilon);
count = count +1;

%25% //30. Fisher equation
yf(count,j) = R*infl - i;
count = count +1;

%26% //31. Interest rate rule
yf(count,j) = i^rho_i*((1/betta)*infl^kappa_pi*(X/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i) - i;
count = count +1;

j=j+1;
end;
