% Author: Max Lampe
% Last Change: 24.03.2015
% to Do: 
% - function for csolve steady state
% - trnsformed variables back into levels
% - reduce system by extra variables: Keff, w, VMPK, prem, Welf  

function yf=GK_ML_stst_csolve(n)
[rows,cols] = size(n);
j=1;
while j<=cols;
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

nuEZ = 0.33;


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
% varrho = n(count,j); count = count+1;
x = n(count,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=1;
% Stochastic discount rate
Lambda = 1;

% Markup
X=1/Pm;
% Euler Eequation for Capital
R = 1 /(Lambda*betta);
a=1;
ksi=1;
g=1;
z=x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM:


%//1.  Marginal utility of consumption
% yf(count,j) = nuEZ * C^(nuEZ-1) * (1-L)^(1-nuEZ) - varrho;
% count = count +1;


%//2.  Marginal utility of labor
% yf(count,j) = (1-nuEZ) * C^(nuEZ) * (1-L)^(-nuEZ) + varrho * Pm*(1-alfa)*Y/L;
% count = count +1;


%//2. Euler equation
% yf(count,j) = betta*exp(R)*exp(Lambda(+1)) -1;

%//3. Stochastic discount rate
% yf(count,j) =   exp(varrho)/exp(varrho(-1)) -  exp(Lambda) ;

yf(count,j) = C/(1-L) *(1-nuEZ)/nuEZ - Pm*(1-alfa)*Y/L;
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
