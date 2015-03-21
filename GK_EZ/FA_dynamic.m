function [residual, g1, g2, g3] = FA_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           columns: equations in order of declaration
%                                                           rows: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              columns: equations in order of declaration
%                                                              rows: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(34, 1);
T20 = exp(y(18))^params(34)*(1-exp(y(16)))^(1-params(34));
T24 = (1-params(36))/params(37);
T32 = (1-params(1))*T20^T24+params(1)*(y(60)^(1-params(36)))^(1/params(37));
T38 = exp(y(18))*(1-params(34))/params(34)/(1-exp(y(16)));
T58 = exp(y(48))^params(34)*(1-exp(y(47)))^(1-params(34));
T60 = (T58/T20)^T24;
T61 = exp(y(18))/exp(y(48));
T62 = T60*T61;
T92 = params(1)*(1-params(5))*exp(y(49))*(exp(y(50))-exp(y(23)))+params(5)*params(1)*exp(y(49))*exp(y(55))*exp(y(52));
T169 = exp(y(32))*params(6)*exp(y(14))/exp(y(1))+exp(y(42))*(exp(y(20))-exp(y(39)));
T181 = exp(y(41))*(exp(y(1))*exp(y(42))*exp(y(33)))^params(6);
T193 = (y(40)+params(33))/(params(33)+y(8))-1;
T195 = params(8)/2*T193^2;
T206 = params(1)*exp(y(49))*params(8)*((params(33)+y(58))/(y(40)+params(33))-1);
T207 = ((params(33)+y(58))/(y(40)+params(33)))^2;
T208 = T206*T207;
T220 = exp(y(32))*params(6)*exp(y(14))/exp(y(33));
T224 = exp(y(1))*exp(y(42))*params(30)*exp(y(33))^params(4);
T262 = params(10)*exp(y(6))*exp(y(12))^((-params(11))*params(9));
T266 = T262*exp(y(44))^params(9);
T275 = (1-params(10)*exp(y(12))^(params(11)*(1-params(10)))*exp(y(44))^(params(10)-1))/(1-params(10));
T300 = exp(y(49))*params(1)*params(10)*exp(y(59))^params(9)*exp(y(44))^(params(11)*(-params(9)))*exp(y(56));
T307 = exp(y(49))*params(1)*params(10)*exp(y(59))^(params(9)-1);
T314 = T307*exp(y(44))^(params(11)*(1-params(9)))*exp(y(57));
T322 = exp(y(44))*exp(y(36))*params(9)/(params(9)-1)/exp(y(37));
T338 = exp(y(7))^params(14);
T342 = 1/params(1)*exp(y(44))^params(12);
T343 = exp(y(34))/(params(9)/(params(9)-1));
T348 = (T342*T343^params(13))^(1-params(14));
T372 = (-(exp(y(32))*(1-params(6))*exp(y(13))/exp(y(16))));
T377 = (-(exp(y(32))*params(6)*exp(y(14))/exp(y(1))/exp(y(2))));
T389 = (-(exp(y(16))^(1-params(6))*exp(y(41))*exp(y(1))*exp(y(42))*exp(y(33))*getPowerDeriv(exp(y(1))*exp(y(42))*exp(y(33)),params(6),1)));
T396 = exp(y(18))^params(34)*(-exp(y(16)))*getPowerDeriv(1-exp(y(16)),1-params(34),1);
T397 = getPowerDeriv(T20,T24,1);
T400 = getPowerDeriv(T32,params(37)/(1-params(36)),1);
T416 = getPowerDeriv(T58/T20,T24,1);
T435 = (1-exp(y(16)))^(1-params(34))*exp(y(18))*getPowerDeriv(exp(y(18)),params(34),1);
T528 = getPowerDeriv(T342*T343^params(13),1-params(14),1);
T555 = params(8)/2*(-(y(40)+params(33)))/((params(33)+y(8))*(params(33)+y(8)))*2*T193;
T567 = params(8)/2*2*T193*1/(params(33)+y(8));
T611 = getPowerDeriv(T275,(-params(9))/(1-params(10)),1);
lhs =y(46);
rhs =T32^(params(37)/(1-params(36)));
residual(1)= lhs-rhs;
lhs =T38;
rhs =exp(y(32))*(1-params(6))*exp(y(13))/exp(y(16));
residual(2)= lhs-rhs;
lhs =exp(y(21));
rhs =T62;
residual(3)= lhs-rhs;
lhs =params(1)*exp(y(21))*exp(y(51));
rhs =1;
residual(4)= lhs-rhs;
lhs =exp(y(27));
rhs =T92;
residual(5)= lhs-rhs;
lhs =exp(y(28));
rhs =1-params(5)+params(5)*params(1)*exp(y(49))*exp(y(54))*exp(y(53));
residual(6)= lhs-rhs;
lhs =exp(y(29));
rhs =exp(y(28))/(params(28)-exp(y(27)));
residual(7)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(3))+(exp(y(22))-exp(y(3)))*exp(y(5));
residual(8)= lhs-rhs;
lhs =exp(y(31));
rhs =exp(y(30))*exp(y(29))/exp(y(5));
residual(9)= lhs-rhs;
lhs =exp(y(20))*exp(y(15));
rhs =exp(y(29))*exp(y(24));
residual(10)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(25))+exp(y(26));
residual(11)= lhs-rhs;
lhs =exp(y(25));
rhs =params(5)*exp(y(30))*exp(y(4))*exp((-x(it_, 4)));
residual(12)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(20))*params(27)*exp(y(42))*exp(y(1));
residual(13)= lhs-rhs;
lhs =exp(y(22));
rhs =T169/exp(y(2));
residual(14)= lhs-rhs;
lhs =exp(y(14));
rhs =T181*exp(y(16))^(1-params(6));
residual(15)= lhs-rhs;
lhs =exp(y(20));
rhs =1+T195+(y(40)+params(33))*params(8)*T193/(params(33)+y(8))-T208;
residual(16)= lhs-rhs;
lhs =exp(y(39));
rhs =params(31)+params(30)/(1+params(4))*exp(y(33))^(1+params(4));
residual(17)= lhs-rhs;
lhs =T220;
rhs =T224;
residual(18)= lhs-rhs;
lhs =y(40);
rhs =exp(y(17))-exp(y(1))*exp(y(42))*exp(y(39));
residual(19)= lhs-rhs;
lhs =exp(y(15));
rhs =y(40)+exp(y(42))*exp(y(1));
residual(20)= lhs-rhs;
lhs =exp(y(19));
rhs =params(32)*exp(y(43));
residual(21)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(17))+exp(y(18))+exp(y(19))+(y(40)+params(33))*T195;
residual(22)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(13))*exp(y(35));
residual(23)= lhs-rhs;
lhs =exp(y(35));
rhs =T266+(1-params(10))*T275^((-params(9))/(1-params(10)));
residual(24)= lhs-rhs;
lhs =exp(y(34));
rhs =1/exp(y(32));
residual(25)= lhs-rhs;
lhs =exp(y(36));
rhs =exp(y(32))*exp(y(13))+T300;
residual(26)= lhs-rhs;
lhs =exp(y(37));
rhs =exp(y(13))+T314;
residual(27)= lhs-rhs;
lhs =exp(y(45));
rhs =T322;
residual(28)= lhs-rhs;
lhs =exp(y(44))^(1-params(9));
rhs =params(10)*exp(y(12))^(params(11)*(1-params(9)))+(1-params(10))*exp(y(45))^(1-params(9));
residual(29)= lhs-rhs;
lhs =exp(y(38));
rhs =exp(y(23))*exp(y(59));
residual(30)= lhs-rhs;
lhs =exp(y(38));
rhs =T338*T348*exp(x(it_, 5));
residual(31)= lhs-rhs;
lhs =y(41);
rhs =params(17)*y(9)-x(it_, 1);
residual(32)= lhs-rhs;
lhs =y(42);
rhs =params(15)*y(10)-x(it_, 2);
residual(33)= lhs-rhs;
lhs =y(43);
rhs =params(19)*y(11)-x(it_, 3);
residual(34)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(34, 65);

  %
  % Jacobian matrix
  %

  g1(1,16)=(-((1-params(1))*T396*T397*T400));
  g1(1,18)=(-(T400*(1-params(1))*T397*T435));
  g1(1,46)=1;
  g1(1,60)=(-(T400*params(1)*getPowerDeriv(y(60),1-params(36),1)*getPowerDeriv(y(60)^(1-params(36)),1/params(37),1)));
  g1(2,13)=T372;
  g1(2,16)=(-(exp(y(18))*(1-params(34))/params(34)*(-exp(y(16)))))/((1-exp(y(16)))*(1-exp(y(16))))-(-(exp(y(16))*exp(y(32))*(1-params(6))*exp(y(13))))/(exp(y(16))*exp(y(16)));
  g1(2,18)=T38;
  g1(2,32)=T372;
  g1(3,16)=(-(T61*(-(T58*T396))/(T20*T20)*T416));
  g1(3,47)=(-(T61*T416*exp(y(48))^params(34)*(-exp(y(47)))*getPowerDeriv(1-exp(y(47)),1-params(34),1)/T20));
  g1(3,18)=(-(T62+T61*T416*(-(T58*T435))/(T20*T20)));
  g1(3,48)=(-(T61*T416*(1-exp(y(47)))^(1-params(34))*exp(y(48))*getPowerDeriv(exp(y(48)),params(34),1)/T20+T60*(-(exp(y(18))*exp(y(48))))/(exp(y(48))*exp(y(48)))));
  g1(3,21)=exp(y(21));
  g1(4,21)=params(1)*exp(y(21))*exp(y(51));
  g1(4,51)=params(1)*exp(y(21))*exp(y(51));
  g1(5,49)=(-T92);
  g1(5,50)=(-(params(1)*(1-params(5))*exp(y(49))*exp(y(50))));
  g1(5,23)=(-(params(1)*(1-params(5))*exp(y(49))*(-exp(y(23)))));
  g1(5,27)=exp(y(27));
  g1(5,52)=(-(params(5)*params(1)*exp(y(49))*exp(y(55))*exp(y(52))));
  g1(5,55)=(-(params(5)*params(1)*exp(y(49))*exp(y(55))*exp(y(52))));
  g1(6,49)=(-(params(5)*params(1)*exp(y(49))*exp(y(54))*exp(y(53))));
  g1(6,28)=exp(y(28));
  g1(6,53)=(-(params(5)*params(1)*exp(y(49))*exp(y(54))*exp(y(53))));
  g1(6,54)=(-(params(5)*params(1)*exp(y(49))*exp(y(54))*exp(y(53))));
  g1(7,27)=(-((-(exp(y(28))*(-exp(y(27)))))/((params(28)-exp(y(27)))*(params(28)-exp(y(27))))));
  g1(7,28)=(-(exp(y(28))/(params(28)-exp(y(27)))));
  g1(7,29)=exp(y(29));
  g1(8,22)=(-(exp(y(22))*exp(y(5))));
  g1(8,3)=(-(exp(y(3))+exp(y(5))*(-exp(y(3)))));
  g1(8,5)=(-((exp(y(22))-exp(y(3)))*exp(y(5))));
  g1(8,30)=exp(y(30));
  g1(9,5)=(-(exp(y(30))*(-(exp(y(29))*exp(y(5))))/(exp(y(5))*exp(y(5)))));
  g1(9,29)=(-(exp(y(30))*exp(y(29))/exp(y(5))));
  g1(9,30)=(-(exp(y(30))*exp(y(29))/exp(y(5))));
  g1(9,31)=exp(y(31));
  g1(10,15)=exp(y(20))*exp(y(15));
  g1(10,20)=exp(y(20))*exp(y(15));
  g1(10,24)=(-(exp(y(29))*exp(y(24))));
  g1(10,29)=(-(exp(y(29))*exp(y(24))));
  g1(11,24)=exp(y(24));
  g1(11,25)=(-exp(y(25)));
  g1(11,26)=(-exp(y(26)));
  g1(12,4)=(-(params(5)*exp(y(30))*exp(y(4))*exp((-x(it_, 4)))));
  g1(12,25)=exp(y(25));
  g1(12,30)=(-(params(5)*exp(y(30))*exp(y(4))*exp((-x(it_, 4)))));
  g1(12,64)=(-(params(5)*exp(y(30))*exp(y(4))*(-exp((-x(it_, 4))))));
  g1(13,1)=(-(exp(y(20))*params(27)*exp(y(42))*exp(y(1))));
  g1(13,20)=(-(exp(y(20))*params(27)*exp(y(42))*exp(y(1))));
  g1(13,26)=exp(y(26));
  g1(13,42)=(-(exp(y(20))*params(27)*exp(y(42))*exp(y(1))));
  g1(14,14)=T377;
  g1(14,1)=(-((-(exp(y(1))*exp(y(32))*params(6)*exp(y(14))))/(exp(y(1))*exp(y(1)))/exp(y(2))));
  g1(14,2)=(-((-(T169*exp(y(2))))/(exp(y(2))*exp(y(2)))));
  g1(14,20)=(-(exp(y(20))*exp(y(42))/exp(y(2))));
  g1(14,22)=exp(y(22));
  g1(14,32)=T377;
  g1(14,39)=(-(exp(y(42))*(-exp(y(39)))/exp(y(2))));
  g1(14,42)=(-(exp(y(42))*(exp(y(20))-exp(y(39)))/exp(y(2))));
  g1(15,14)=exp(y(14));
  g1(15,1)=T389;
  g1(15,16)=(-(T181*exp(y(16))*getPowerDeriv(exp(y(16)),1-params(6),1)));
  g1(15,33)=T389;
  g1(15,41)=(-(T181*exp(y(16))^(1-params(6))));
  g1(15,42)=T389;
  g1(16,20)=exp(y(20));
  g1(16,49)=T208;
  g1(16,8)=(-(T555+((params(33)+y(8))*(y(40)+params(33))*params(8)*(-(y(40)+params(33)))/((params(33)+y(8))*(params(33)+y(8)))-(y(40)+params(33))*params(8)*T193)/((params(33)+y(8))*(params(33)+y(8)))));
  g1(16,40)=(-(T567+(params(8)*T193+(y(40)+params(33))*params(8)*1/(params(33)+y(8)))/(params(33)+y(8))-(T207*params(1)*exp(y(49))*params(8)*(-(params(33)+y(58)))/((y(40)+params(33))*(y(40)+params(33)))+T206*(-(params(33)+y(58)))/((y(40)+params(33))*(y(40)+params(33)))*2*(params(33)+y(58))/(y(40)+params(33)))));
  g1(16,58)=T207*params(1)*exp(y(49))*params(8)*1/(y(40)+params(33))+T206*2*(params(33)+y(58))/(y(40)+params(33))*1/(y(40)+params(33));
  g1(17,33)=(-(params(30)/(1+params(4))*exp(y(33))*getPowerDeriv(exp(y(33)),1+params(4),1)));
  g1(17,39)=exp(y(39));
  g1(18,14)=T220;
  g1(18,1)=(-T224);
  g1(18,32)=T220;
  g1(18,33)=(-(exp(y(32))*params(6)*exp(y(14))*exp(y(33))))/(exp(y(33))*exp(y(33)))-exp(y(1))*exp(y(42))*params(30)*exp(y(33))*getPowerDeriv(exp(y(33)),params(4),1);
  g1(18,42)=(-T224);
  g1(19,1)=exp(y(1))*exp(y(42))*exp(y(39));
  g1(19,17)=(-exp(y(17)));
  g1(19,39)=exp(y(1))*exp(y(42))*exp(y(39));
  g1(19,40)=1;
  g1(19,42)=exp(y(1))*exp(y(42))*exp(y(39));
  g1(20,1)=(-(exp(y(42))*exp(y(1))));
  g1(20,15)=exp(y(15));
  g1(20,40)=(-1);
  g1(20,42)=(-(exp(y(42))*exp(y(1))));
  g1(21,19)=exp(y(19));
  g1(21,43)=(-(params(32)*exp(y(43))));
  g1(22,13)=exp(y(13));
  g1(22,17)=(-exp(y(17)));
  g1(22,18)=(-exp(y(18)));
  g1(22,19)=(-exp(y(19)));
  g1(22,8)=(-((y(40)+params(33))*T555));
  g1(22,40)=(-(T195+(y(40)+params(33))*T567));
  g1(23,13)=(-(exp(y(13))*exp(y(35))));
  g1(23,14)=exp(y(14));
  g1(23,35)=(-(exp(y(13))*exp(y(35))));
  g1(24,6)=(-T266);
  g1(24,35)=exp(y(35));
  g1(24,12)=(-(exp(y(44))^params(9)*params(10)*exp(y(6))*exp(y(12))*getPowerDeriv(exp(y(12)),(-params(11))*params(9),1)+(1-params(10))*(-(exp(y(44))^(params(10)-1)*params(10)*exp(y(12))*getPowerDeriv(exp(y(12)),params(11)*(1-params(10)),1)))/(1-params(10))*T611));
  g1(24,44)=(-(T262*exp(y(44))*getPowerDeriv(exp(y(44)),params(9),1)+(1-params(10))*T611*(-(params(10)*exp(y(12))^(params(11)*(1-params(10)))*exp(y(44))*getPowerDeriv(exp(y(44)),params(10)-1,1)))/(1-params(10))));
  g1(25,32)=(-((-exp(y(32)))/(exp(y(32))*exp(y(32)))));
  g1(25,34)=exp(y(34));
  g1(26,13)=(-(exp(y(32))*exp(y(13))));
  g1(26,49)=(-T300);
  g1(26,32)=(-(exp(y(32))*exp(y(13))));
  g1(26,36)=exp(y(36));
  g1(26,56)=(-T300);
  g1(26,44)=(-(exp(y(56))*exp(y(49))*params(1)*params(10)*exp(y(59))^params(9)*exp(y(44))*getPowerDeriv(exp(y(44)),params(11)*(-params(9)),1)));
  g1(26,59)=(-(exp(y(56))*exp(y(44))^(params(11)*(-params(9)))*exp(y(49))*params(1)*params(10)*exp(y(59))*getPowerDeriv(exp(y(59)),params(9),1)));
  g1(27,13)=(-exp(y(13)));
  g1(27,49)=(-T314);
  g1(27,37)=exp(y(37));
  g1(27,57)=(-T314);
  g1(27,44)=(-(exp(y(57))*T307*exp(y(44))*getPowerDeriv(exp(y(44)),params(11)*(1-params(9)),1)));
  g1(27,59)=(-(exp(y(57))*exp(y(44))^(params(11)*(1-params(9)))*exp(y(49))*params(1)*params(10)*exp(y(59))*getPowerDeriv(exp(y(59)),params(9)-1,1)));
  g1(28,36)=(-T322);
  g1(28,37)=(-(exp(y(44))*(-(exp(y(37))*exp(y(36))*params(9)/(params(9)-1)))/(exp(y(37))*exp(y(37)))));
  g1(28,44)=(-T322);
  g1(28,45)=exp(y(45));
  g1(29,12)=(-(params(10)*exp(y(12))*getPowerDeriv(exp(y(12)),params(11)*(1-params(9)),1)));
  g1(29,44)=exp(y(44))*getPowerDeriv(exp(y(44)),1-params(9),1);
  g1(29,45)=(-((1-params(10))*exp(y(45))*getPowerDeriv(exp(y(45)),1-params(9),1)));
  g1(30,23)=(-(exp(y(23))*exp(y(59))));
  g1(30,38)=exp(y(38));
  g1(30,59)=(-(exp(y(23))*exp(y(59))));
  g1(31,34)=(-(exp(x(it_, 5))*T338*T342*T343*getPowerDeriv(T343,params(13),1)*T528));
  g1(31,7)=(-(exp(x(it_, 5))*T348*exp(y(7))*getPowerDeriv(exp(y(7)),params(14),1)));
  g1(31,38)=exp(y(38));
  g1(31,44)=(-(exp(x(it_, 5))*T338*T528*T343^params(13)*1/params(1)*exp(y(44))*getPowerDeriv(exp(y(44)),params(12),1)));
  g1(31,65)=(-(T338*T348*exp(x(it_, 5))));
  g1(32,9)=(-params(17));
  g1(32,41)=1;
  g1(32,61)=1;
  g1(33,10)=(-params(15));
  g1(33,42)=1;
  g1(33,62)=1;
  g1(34,11)=(-params(19));
  g1(34,43)=1;
  g1(34,63)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],34,4225);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],34,274625);
end
end
