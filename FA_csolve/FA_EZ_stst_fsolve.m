% Author: Max Lampe
% Last Change: 07.03.2015
% to Do: 
% - function for csolve steady state
% - trnsformed variables back into levels
% - reduce system by extra variables: Keff, w, VMPK, prem, Welf  

function yf=FA_EZ_stst_fsolve(n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% parameters in the system:

betta=0.99000000;
sig=1.00000000;
varphi=0.27600000;
zetta=7.20000000;
theta=0.97155955;
alfa=0.33000000;
G_over_Y=0.20000000;
eta_i=1.72800000;
epsilon=4.16700000;
gam=0.77900000;
gam_P=0.24100000;
kappa_pi=1.50000000;
kappa_y=-0.12500000;
rho_i=0.00000000;
rho_ksi=0.66000000;
sigma_ksi=0.05000000;
rho_a=0.95000000;
sigma_a=0.01000000;
rho_g=0.95000000;
sigma_g=0.01000000;
sigma_Ne=0.01000000;
sigma_i=0.01000000;
rho_shock_psi=0.66000000;
sigma_psi=0.07200000;
kappa=10.00000000;
tau=0.00100000;
omega=0.00222778;
lambda=0.38149499;
chi=3.41080850;
b=0.03760101;
delta_c=0.02041451;
G_ss       =   0.16975710;
I_ss       =   0.14153927;

gammaEZ = 2;
psiEZ = 0.5;
thetaEZ = (1-gammaEZ)/(1-1/psiEZ);
nuEZ = 0.8;


%%%%%%%%%%%%% reducing the system of equation %%%%%%%%%%%%% 


% variables in the system:
count=1;
C = n(count); count = count+1;
D = n(count); count = count+1;
F = n(count); count = count+1;
G = n(count); count = count+1;
I = n(count); count = count+1;
In = n(count); count = count+1;
K = n(count); count = count+1;
L = n(count); count = count+1;
N = n(count); count = count+1;
Ne = n(count); count = count+1;
Nn = n(count); count = count+1;
Pm = n(count); count = count+1;
Q = n(count); count = count+1;
Rk = n(count); count = count+1;
U = n(count); count = count+1;
Y = n(count); count = count+1;
Ym = n(count); count = count+1;
Z = n(count); count = count+1;
delta = n(count); count = count+1;
eta = n(count); count = count+1;
i = n(count); count = count+1;
infl = n(count); count = count+1;
inflstar = n(count); count = count+1;
nu = n(count); count = count+1;
phi = n(count); count = count+1;
v = n(count); count = count+1;
x = n(count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
count=1;
% EZ Zin Preferences
%1% // 0.3 Value Function
yf(count) = ((1-betta)*u^((1-gammaEZ)/thetaEZ) + betta*((v^(1-gammaEZ))^(1/thetaEZ)))^(thetaEZ/(1-gammaEZ))- v;
count = count +1;

%2% // 0.4 Static Leisure-Consumption
yf(count) = Pm*(1-alfa)*Y/L - (1-nuEZ)/nuEZ*C/(1-L);
count = count +1;

% //Financial Intermediaries
%5%  //5. Value of banks' capital
 yf(count) = (1-theta)*betta*Lambda*(Rk-R)+betta*Lambda*theta*x*nu-nu;
count = count +1;

% //6. Value of banks' net wealth
yf(count) = (1-theta) + betta * Lambda * theta * z * eta - eta;
count = count +1;

% //7. Optimal leverage
yf(count) = eta/(lambda-nu) - phi;
count = count +1;

%6% //8. Growth rate of banks' capital
yf(count) = (Rk-R) * phi + R - z;
count = count +1;

% //9. Growth rate of banks' net wealth
% yf(count) = z-x;
% count = count +1;

% //Aggregate capital, net worth
%7% //10. Aggregate capital
yf(count) = phi * N - Q * K;
count = count +1;

%8% //11. Banks' net worth
yf(count) = Ne+Nn-N;
count = count +1;

%9% //12. Existing banks' net worth accumulation
yf(count) = theta*z*N-Ne;
count = count +1;

% //13. New banks' net worth
yf(count) = omega*Q*ksi*K-Nn;
count = count +1;

% //Final goods producer
%10% //14. Return to capital
yf(count) = (Pm * alfa * Ym/K + ksi * (Q-delta) )/ Q - Rk;
count = count +1;

%11% //15. Production function
yf(count) = a*(ksi*U*K)^alfa*L^(1-alfa) - Ym;
count = count +1;

% //Capital Goods Producer
%12% //16. Optimal investment decision
yf(count) = 1+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2 ...
+eta_i*((In+I_ss)/(In+I_ss)-1)*(In+I_ss)/(In+I_ss)-betta*Lambda*eta_i...
*((In+I_ss)/(In+I_ss)-1)*((In+I_ss)/(In+I_ss))^2 - Q;
count = count +1;

%13% //17. Depreciation rate
yf(count) = delta_c+b/(1+zetta)*U^(1+zetta) - delta;
count = count +1;

%14% //18. Optimal capacity utilization rate
yf(count) = b*U^zetta*ksi*K - Pm*alfa*Ym/U;
count = count +1;

%15% //19. Net investment
yf(count) = I-delta*ksi*K - In;
count = count +1;

%16% //20. Capital accumulation equation
yf(count) = ksi*K+In-K; 
count = count +1;

%17% //21. Government consumption
yf(count) = G_ss*g-G;
count = count +1;

% //Equilibrium
%18% //22. Aggregate resource constraint
yf(count) = C+G+I+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2*(In+I_ss) - Y;
count = count +1;

%19% //23. Wholesale, retail output
yf(count) = Y*D-Ym;
count = count +1;

%20% //24. Price dispersion
yf(count) = gam*D*infl^(-gam_P*epsilon)*infl^epsilon ...
+(1-gam)*((1-gam*infl^(gam_P*(1-gam))*infl^(gam-1))/(1-gam))^(-epsilon/(1-gam)) - D;
count = count +1;

% //25. Markup
% yf(count) = 1/Pm-X;
% count = count +1;

%21% //26. Optimal price choice
yf(count) = Y*Pm+betta*gam*Lambda*infl^epsilon*(infl)^(-epsilon*gam_P)*F- F;
count = count +1;

%22% //27.
yf(count) = Y+betta*gam*Lambda*infl^(epsilon-1)*infl^(gam_P*(1-epsilon))*Z-Z;
count = count +1;

%23% //28. Optimal price choice
yf(count) = epsilon/(epsilon-1)*F/Z*infl-inflstar;
count = count +1;

%24% //29. Price index
yf(count) = gam*infl^(gam_P*(1-epsilon))+(1-gam)*(inflstar)^(1-epsilon) - infl^(1-epsilon);
count = count +1;

%25% //30. Fisher equation
yf(count) = R*infl - i;
count = count +1;

%26% //31. Interest rate rule
yf(count) = i^rho_i*((1/betta)*infl^kappa_pi*(X/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i) - i;




