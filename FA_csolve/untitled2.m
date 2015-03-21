
count = 1;
Y = n(count,j); count = count+1;            %1 
Ym = n(count,j); count = count+1;            %2
K = n(count,j); count = count+1;            %3
% Keff = n(4,j);
L = n(count,j); count = count+1;            %4
I = n(count,j); count = count+1;            %5
C = n(count,j); count = count+1;            %6
G = n(count,j); count = count+1;            %7
Q = n(count,j); count = count+1;            %8
varrho = n(count,j); count = count+1;            %9
Lambda = n(count,j); count = count+1;            %10
Rk = n(count,j); count = count+1;            %11
R = n(count,j); count = count+1;            %12
N = n(count,j); count = count+1;            %13
Ne = n(count,j); count = count+1;            %14
Nn = n(count,j); count = count+1;            %15
nu = n(count,j); count = count+1;            %16
eta = n(count,j); count = count+1;            %17
phi = n(count,j); count = count+1;            %18
z = n(count,j); count = count+1;            %20
x = n(count,j); count = count+1;            %21
Pm = n(count,j); count = count+1;            %22
% w = n(23,j); count 
% VMPK = n(24,j);
U = n(count,j); count = count+1;            %23
X = n(count,j); count = count+1;            %24
D = n(count,j); count = count+1;            %25
F = n(count,j); count = count+1;            %26
Z = n(count,j); count = count+1;            %27
i = n(count,j); count = count+1;            %28
% prem = n(31,j);
delta = n(32,j); count = count+1;           %29
In = n(count,j); count = count+1;           %30
% Welf = n(34,j);
a = n(count,j); count = count+1;            %31
ksi = n(count,j); count = count+1;          %32
g = n(count,j); count = count+1;            %33
infl = n(count,j); count = count+1;         %34
inflstar = n(count,j);                      %35



%% level different initial Values for Steady State Optimization without extra variables
initial = [0.9*ones(32,1); zeros(2,1); 0.9; 0.9; 0.9 ;0.9;];
initial = [0.9*ones(29,1); zeros(2,1); 0.9; 0.9; 0.9 ;0.9;];
% initial = 0.9*ones(38,1)
initial = 0.9*ones(35,1)

[SS_opti, rc]=csolve(@FA_stst_csolve_2_trans,initial,[],1e-12,3)

% mean(abs(SS_opti - m'))

% [SS_opti, rc]=csolve(@FA_stst_csolve_exp_2 ,initial,[],1e-12,10000);
% mean(abs(SS_opti - m'))




%% exponential different initial Values for Steady State Optimization
initial = [0.9*ones(28,1); m(33:end)'];
[SS_opti, rc]=csolve(@FA_stst_csolve_exp_2,initial,[],1e-12,1000);
mean(abs(SS_opti - m'))

[SS_opti, rc]=csolve(@FA_stst_csolve_exp_2 ,initial,[],1e-12,10000);
mean(abs(SS_opti - m'))




%% exponential different initial Values for Steady State Optimization
initial = [0.9*ones(32,1); m(33:end)'];
[SS_opti, rc]=csolve(@FA_stst_csolve_exp,initial,[],1e-12,1000);
mean(abs(SS_opti - m'))

[SS_opti, rc]=csolve(@FA_stst_csolve_exp,initial,[],1e-12,10000);
mean(abs(SS_opti - m'))


%% level different initial Values for Steady State Optimization
initial = [-0.9*ones(32,1); m(33:end)'];
SS_Paper = [exp(m(1:32)'); m(33:end)'];
[SS_opti ,rc]=csolve(@FA_stst_csolve_2,initial,[],1e-12,10);
mean(abs(SS_opti - SS_Paper))

[SS_opti ,rc]=csolve(@FA_stst_csolve_2,initial,[],1e-12,10000);
mean(abs(SS_opti - SS_Paper))

 
%% levels steady as initial
initial = exp(values);
initial(end-5) =  -296.6207;
[SS_VALUES_paper,rc]=csolve(@FA_stst_csolve_2,initial,[],1e-12,1000)
% max(SS_VALUES_paper - x1)


