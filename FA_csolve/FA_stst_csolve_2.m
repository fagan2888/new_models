function yf=FA_stst_csolve_2(n)
[rows,cols] = size(n);
j=1;
while j<=cols;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% parameters in the system:

betta=0.99000000;
sig=1.00000000;
hh=0.81500000;
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

% variables in the system:
Y = n(1,j);
Ym = n(2,j);
K = n(3,j);
Keff = n(4,j);
L = n(5,j);
I = n(6,j);
C = n(7,j);
G = n(8,j);
Q = n(9,j);
varrho = n(10,j);
Lambda = n(11,j);
Rk = n(12,j);
R = n(13,j);
N = n(14,j);
Ne = n(15,j);
Nn = n(16,j);
nu = n(17,j);
eta = n(18,j);
phi = n(19,j);
z = n(20,j);
x = n(21,j);
Pm = n(22,j);
w = n(23,j);
VMPK = n(24,j);
U = n(25,j);
X = n(26,j);
D = n(27,j);
F = n(28,j);
Z = n(29,j);
i = n(30,j);
prem = n(31,j);
delta = n(32,j);
In = n(33,j);
Welf = n(34,j);
a = n(35,j);
ksi = n(36,j);
g = n(37,j);
infl = n(38,j);
inflstar = n(39,j);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM:
% //1. Marginal utility of consumption
yf(1,j) = (C-hh*C)^(-sig)-betta*hh*(C-hh*C)^(-sig) - varrho;

% //2. Euler equation
yf(2,j) = 1 - betta*R*Lambda;

% //3. Stochastic discount rate
yf(3,j) = varrho/varrho - Lambda;

% //4. Labor market equilibrium
yf(4,j) = varrho*Pm*(1-alfa)*Y/L - chi*L^varphi;

% //Financial Intermediaries
% //5. Value of banks' capital
 yf(5,j) = (1-theta)*betta*Lambda*(Rk-R)+betta*Lambda*theta*x*nu-nu;

% //6. Value of banks' net wealth
yf(6,j) = (1-theta)+betta*Lambda*theta*z*eta-eta;

% //7. Optimal leverage
yf(7,j) = eta/(lambda-nu)- phi;

% //8. Growth rate of banks' capital
yf(8,j) = (Rk-R)*phi+R-z;

% //9. Growth rate of banks' net wealth
yf(9,j) = (phi/phi)*z-x;

% //Aggregate capital, net worth
% //10. Aggregate capital
yf(10,j) = phi*N-Q*K;

% //11. Banks' net worth
yf(11,j) = Ne+Nn-N;

% //12. Existing banks' net worth accumulation
yf(12,j) = theta*z*N-Ne;

% //13. New banks' net worth
yf(13,j) = omega*Q*ksi*K-Nn;

% //Final goods producer
% //14. Return to capital
yf(14,j) = (Pm*alfa*Ym/K+ksi*(Q-delta))/Q-Rk;

% //15. Production function
yf(15,j) = a*(ksi*U*K)^alfa*L^(1-alfa)-Ym;

% //Capital Goods Producer
% //16. Optimal investment decision
yf(16,j) = 1+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2 ...
+eta_i*((In+I_ss)/(In+I_ss)-1)*(In+I_ss)/(In+I_ss)-betta*Lambda*eta_i...
*((In+I_ss)/(In+I_ss)-1)*((In+I_ss)/(In+I_ss))^2 - Q;

% //17. Depreciation rate
yf(17,j) = delta_c+b/(1+zetta)*U^(1+zetta) - delta;

% //18. Optimal capacity utilization rate
yf(18,j) = b*U^zetta*ksi*K - Pm*alfa*Ym/U;

% //19. Net investment
yf(19,j) = I-delta*ksi*K - In;

% //20. Capital accumulation equation
yf(20,j) = ksi*K+In-K; 

% //21. Government consumption
yf(21,j) = G_ss*g-G;

% //Equilibrium
% //22. Aggregate resource constraint
yf(22,j) = C+G+I+eta_i/2*((In+I_ss)/(In+I_ss)-1)^2*(In+I_ss) - Y;

% //23. Wholesale, retail output
yf(23,j) = Y*D-Ym;

% //24. Price dispersion
yf(24,j) = gam*D*infl^(-gam_P*epsilon)*infl^epsilon ...
+(1-gam)*((1-gam*infl^(gam_P*(1-gam))*infl^(gam-1))/(1-gam))^(-epsilon/(1-gam)) - D;

% //25. Markup
yf(25,j) = 1/Pm-X;

% //26. Optimal price choice
yf(26,j) = Y*Pm+betta*gam*Lambda*infl^epsilon*(infl)^(-epsilon*gam_P)*F(+1)- F;

% //27.
yf(27,j) = Y+betta*gam*Lambda*infl^(epsilon-1)*infl^(gam_P*(1-epsilon))*Z-Z;

% //28. Optimal price choice
yf(28,j) = epsilon/(epsilon-1)*F/Z*infl-inflstar;

% //29. Price index
yf(29,j) = gam*infl^(gam_P*(1-epsilon))+(1-gam)*(inflstar)^(1-epsilon) - (infl)^(1-epsilon);

% //30. Fisher equation
yf(30,j) = R*infl - i;

% //31. Interest rate rule
yf(31,j) = i^rho_i*((1/betta)*infl^kappa_pi*(X/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i) - i;

% //Shocks
% //32. TFP shock
yf(32,j) = rho_a*a - a;

% //33. Capital quality shock
yf(33,j) = rho_ksi*ksi - ksi;

% //34. Government consumption shock
yf(34,j) = rho_g*g - g;

% //Some extra variables for convenience
% //35. Effective capital
yf(35,j) = ksi*K - Keff;

% //36. Wages
yf(36,j) = Pm*(1-alfa)*Y/L - w;

% //37. Marginal value product of capital
yf(37,j) = Pm*alfa*Y/(ksi*K) - VMPK;


% //38. Welfare
yf(38,j) =   log(C-hh*C)-chi*L^(1+varphi)/(1+varphi)+betta*Welf - Welf;

% //39. Premium
yf(39,j) = Rk/R - prem;

j=j+1;
end;

% call from Matlab Command Window (or from within some m-file by defining a starting vector x0
% 
% x0 = ones(39,1)*0.9;
%
% and calling the fsolve-function on the system defined by yf:
%
% [SS_VALUES,rc]=csolve(@FA_stst_csolve,x0,[],1e-12,1000), or:
% [SS_VALUES,rc]=csolve(@neoclassical_stst_csolve,x0,[],1e-12,1000,BETTA,SIG,ALFA,DELTA) 
% if you pass through the parameters to neoclassical_stst_fsolve.m
%
% type help csolve to find out how to modify the options
