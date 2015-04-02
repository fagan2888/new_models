// Max Lampe
// Last change: 03.003.2015
// Change:
// - exchange utility fucntion Household
// - Epstein Zin rekursive utility function
// - change of Initial values defined in the initial block
// - Blanchard conditions fullfiled
// - 


parameters betta sig varphi zetta theta alfa G_over_Y eta_i epsilon gam gam_P kappa_pi kappa_y rho_i rho_ksi sigma_ksi rho_a sigma_a rho_g sigma_g sigma_Ne sigma_i rho_shock_psi sigma_psi kappa tau omega lambda chi b delta_c G_ss I_ss nuEZ psiEZ gammaEZ thetaEZ cte;
var Y Ym K L I C G Q Lambda Rk R N Ne Nn nu eta phi z x Pm U X D F Z i delta In a ksi g infl inflstar v u ev;
varexo e_a e_ksi e_g e_Ne e_i ;

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
nuEZ = 0.33;
cte = (0.451610879111627^nuEZ*(1-0.292857916980711)^(1-nuEZ))^(-1);

model;
//Home household Epstein Zin Utility
// 0.3 Value Function

// Eppstein Zin Utility function 
u = cte*exp(C)^nuEZ*(1-exp(L))^(1-nuEZ);
ev = v(+1)^(1-gammaEZ);
v = ((1-betta)*u^((1-gammaEZ)/thetaEZ)+ betta*(ev^(1/thetaEZ)))^(thetaEZ/(1-gammaEZ));


// stochastic discount factor (kernel)
exp(Lambda) = (u(+1)/u)^((1-gammaEZ)/thetaEZ)*(exp(C)/exp(C(+1)))*(v(+1)^(1-gammaEZ)/ev)^(1-1/thetaEZ);

// Eulerequation
betta*exp(Lambda)*exp(R(+1)) = 1;

// Marginal rate for labor and consumption
(1-nuEZ)/nuEZ*exp(C)/(1-exp(L)) = exp(Pm)*(1-alfa)*exp(Y)/exp(L);



//Financial Intermediaries
//5. Value of banks' capital
exp(nu)     =   (1-theta)*betta*exp(Lambda(+1))*(exp(Rk(+1))-exp(R))+betta*exp(Lambda(+1))*theta*exp(x(+1))*exp(nu(+1));

//6. Value of banks' net wealth
exp(eta)    =   (1-theta)+betta*exp(Lambda(+1))*theta*exp(z(+1))*exp(eta(+1));

//7. Optimal leverage
exp(phi)    =   exp(eta)/(lambda-exp(nu));

//8. Growth rate of banks' capital
// exp(z)      =   (exp(Rk)-exp(R))*exp(phi)+exp(R)
exp(z)      =   (exp(Rk)-exp(R(-1)))*exp(phi(-1))+exp(R(-1));

//9. Growth rate of banks' net wealth
// exp(x)      =   (exp(phi)/exp(phi))*exp(z)
exp(x)      =   (exp(phi)/exp(phi(-1)))*exp(z);

//Aggregate capital, net worth
//10. Aggregate capital
exp(Q)*exp(K)     =   exp(phi)*exp(N);

//11. Banks' net worth
exp(N)      =   exp(Ne)+exp(Nn);

//12. Existing banks' net worth accumulation
exp(Ne)     =   theta*exp(z)*exp(N(-1))*exp(-e_Ne);

//13. New banks' net worth
exp(Nn)    =   omega*exp(Q)*exp(ksi)*exp(K(-1));

//Final goods producer
//14. Return to capital
exp(Rk)     =   (exp(Pm)*alfa*exp(Ym)/exp(K(-1))+exp(ksi)*(exp(Q)-exp(delta)))/exp(Q(-1));

//15. Production function
exp(Ym)     =   exp(a)*(exp(ksi)*exp(U)*exp(K(-1)))^alfa*exp(L)^(1-alfa);

//Capital Goods Producer
//16. Optimal investment decision
exp(Q)  =   1+eta_i/2*((In+I_ss)/(In(-1)+I_ss)-1)^2+eta_i*((In+I_ss)/(In(-1)+I_ss)-1)*(In+I_ss)/(In(-1)+I_ss)-betta*exp(Lambda(+1))*eta_i*((In(+1)+I_ss)/(In+I_ss)-1)*((In(+1)+I_ss)/(In+I_ss))^2;

//17. Depreciation rate
exp(delta) = delta_c+b/(1+zetta)*exp(U)^(1+zetta);

//18. Optimal capacity utilization rate
exp(Pm)*alfa*exp(Ym)/exp(U) = b*exp(U)^zetta*exp(ksi)*exp(K(-1));

//19. Net investment
In  =   exp(I)-exp(delta)*exp(ksi)*exp(K(-1));

//20. Capital accumulation equation
exp(K)  =   exp(ksi)*exp(K(-1))+In; 

//21. Government consumption
exp(G)   =   G_ss*exp(g);

//Equilibrium
//22. Aggregate resource constraint
exp(Y)   =   exp(C)+exp(G)+exp(I)+eta_i/2*((In+I_ss)/(In(-1)+I_ss)-1)^2*(In+I_ss);

//23. Wholesale, retail output
exp(Ym)    =   exp(Y)*exp(D);

//24. Price dispersion
exp(D)    =   gam*exp(D(-1))*exp(infl(-1))^(-gam_P*epsilon)*exp(infl)^epsilon+(1-gam)*((1-gam*exp(infl(-1))^(gam_P*(1-gam))*exp(infl)^(gam-1))/(1-gam))^(-epsilon/(1-gam));

//25. Markup
exp(X)  =   1/exp(Pm);

//26. Optimal price choice
exp(F)       =   exp(Y)*exp(Pm)+betta*gam*exp(Lambda(+1))*exp(infl(+1))^epsilon*(exp(infl))^(-epsilon*gam_P)*exp(F(+1));

//27.
exp(Z)       =   exp(Y)+betta*gam*exp(Lambda(+1))*exp(infl(+1))^(epsilon-1)*exp(infl)^(gam_P*(1-epsilon))*exp(Z(+1));

//28. Optimal price choice
exp(inflstar)   =  epsilon/(epsilon-1)*exp(F)/exp(Z)*exp(infl);

//29. Price index
(exp(infl))^(1-epsilon)     =   gam*exp(infl(-1))^(gam_P*(1-epsilon))+(1-gam)*(exp(inflstar))^(1-epsilon);

//30. Fisher equation
exp(i)  =   exp(R)*exp(infl(+1));

//31. Interest rate rule
exp(i)      =   exp(i(-1))^rho_i*((1/betta)*exp(infl)^kappa_pi*(exp(X)/(epsilon/(epsilon-1)))^(kappa_y))^(1-rho_i)*exp(e_i);

//Shocks
//32. TFP shock
a  =   rho_a*a(-1)-e_a;

//33. Capital quality shock
ksi=   rho_ksi*ksi(-1)-e_ksi;

//34. Government consumption shock
g  =   rho_g*g(-1)-e_g;

//Some extra variables for convenience
//35. Effective capital
//exp(Keff)   =   exp(ksi)*exp(K(-1));
//36. Wages
// exp(w)      =   exp(Pm)*(1-alfa)*exp(Y)/exp(L);
//37. Marginal value product of capital
// exp(VMPK)   =   exp(Pm)*alfa*exp(Y)/(exp(ksi)*exp(K(-1)));
//38. Welfare
// Welf   =   log(exp(C)-hh*exp(C(-1)))-chi*exp(L)^(1+varphi)/(1+varphi)+betta*Welf(+1);
//39. Premium
// exp(prem)   =   exp(Rk(+1))/exp(R);
end;


initval;
C = log(0.451610879111627);
D = log(1.000000000000000);
F = log(2.477215009286168);
G = log(0.169757100000000);
I = log(0.169757100000000);
In = 0; 
K = log(4.974107642854883);
L = log(0.292857916980711);
N = log(1.243526844272020);
Ne = log(1.232445626747421);
Nn = log(0.011081217524599);
Pm = log(0.760019198464123);
Q = log(1);
Rk = log(1.012601010275084);
U = log(0.999999997412428);
Y = log(0.745720664846253);
Ym = log(0.745720664846254);
Z = log(3.259411096841005);
delta = log(0.024999998927095);
eta = log(1.511020932460962);
i = log(1.010101010101010);
infl = log(1.00000000000000);
inflstar = log(1.00000000000000);
nu = log(0.003739777068239);
phi = log(4.000000213720197);
v =  log(0.609876260046773);
x = log(1.020101011331604);
Lambda = log(1^((1-gammaEZ)/thetaEZ)*1^(1-1/thetaEZ));
X=-Pm;
R = Lambda-log(betta);
a=log(1);
ksi=log(1);
g=log(1);
z=x;
e_a=0.00000000;
e_ksi=0.00000000;
e_g=0.00000000;
e_Ne=0.00000000;
e_i=0.00000000;

u = 1;
v = 1;
ev = 1^(1-gammaEZ);

end;

shocks;
var e_ksi=sigma_ksi^2;
end;

check;

// steady(solve_algo=0); 

//stoch_simul(order=1, periods=2000, irf=40);
 stoch_simul(PRUNING, order=2, periods=2000, irf=40);
// stoch_simul(PRUNING ,order=2, periods=2000, irf=40);

// Saving the impulse responses
// save ../data/FA_1.mat Y_e_a Y_e_ksi Y_e_g Y_e_Ne Y_e_i Ym_e_a Ym_e_ksi Ym_e_g Ym_e_Ne Ym_e_i K_e_a K_e_ksi K_e_g K_e_Ne K_e_i Keff_e_a Keff_e_ksi Keff_e_g Keff_e_Ne Keff_e_i L_e_a L_e_ksi L_e_g L_e_Ne L_e_i I_e_a I_e_ksi I_e_g I_e_Ne I_e_i C_e_a C_e_ksi C_e_g C_e_Ne C_e_i G_e_a G_e_ksi G_e_g G_e_Ne G_e_i Q_e_a Q_e_ksi Q_e_g Q_e_Ne Q_e_i varrho_e_a varrho_e_ksi varrho_e_g varrho_e_Ne varrho_e_i Lambda_e_a Lambda_e_ksi Lambda_e_g Lambda_e_Ne Lambda_e_i Rk_e_a Rk_e_ksi Rk_e_g Rk_e_Ne Rk_e_i R_e_a R_e_ksi R_e_g R_e_Ne R_e_i N_e_a N_e_ksi N_e_g N_e_Ne N_e_i Ne_e_a Ne_e_ksi Ne_e_g Ne_e_Ne Ne_e_i Nn_e_a Nn_e_ksi Nn_e_g Nn_e_Ne Nn_e_i nu_e_a nu_e_ksi nu_e_g nu_e_Ne nu_e_i eta_e_a eta_e_ksi eta_e_g eta_e_Ne eta_e_i phi_e_a phi_e_ksi phi_e_g phi_e_Ne phi_e_i z_e_a z_e_ksi z_e_g z_e_Ne z_e_i x_e_a x_e_ksi x_e_g x_e_Ne x_e_i Pm_e_a Pm_e_ksi Pm_e_g Pm_e_Ne Pm_e_i w_e_a w_e_ksi w_e_g w_e_Ne w_e_i VMPK_e_a VMPK_e_ksi VMPK_e_g VMPK_e_Ne VMPK_e_i U_e_a U_e_ksi U_e_g U_e_Ne U_e_i X_e_a X_e_ksi X_e_g X_e_Ne X_e_i D_e_a D_e_ksi D_e_g D_e_Ne D_e_i F_e_a F_e_ksi F_e_g F_e_Ne F_e_i Z_e_a Z_e_ksi Z_e_g Z_e_Ne Z_e_i i_e_a i_e_ksi i_e_g i_e_Ne i_e_i prem_e_a prem_e_ksi prem_e_g prem_e_Ne prem_e_i delta_e_a delta_e_ksi delta_e_g delta_e_Ne delta_e_i In_e_a In_e_ksi In_e_g In_e_Ne In_e_i Welf_e_a Welf_e_ksi Welf_e_g Welf_e_Ne Welf_e_i a_e_a a_e_ksi a_e_g a_e_Ne a_e_i ksi_e_a ksi_e_ksi ksi_e_g ksi_e_Ne ksi_e_i g_e_a g_e_ksi g_e_g g_e_Ne g_e_i infl_e_a infl_e_ksi infl_e_g infl_e_Ne infl_e_i inflstar_e_a inflstar_e_ksi inflstar_e_g inflstar_e_Ne inflstar_e_i ;
