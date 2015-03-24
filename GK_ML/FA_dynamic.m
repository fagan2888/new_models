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
T18 = exp(y(34))*(1-params(7))*exp(y(14))/exp(y(17));
T26 = exp(y(19))/(1-exp(y(17)))*(1-params(35))/params(35);
T68 = exp(y(48))*params(1)*(1-params(6))*(exp(y(49))-exp(y(25)))+params(6)*params(1)*exp(y(48))*exp(y(53))*exp(y(50));
T145 = exp(y(34))*params(7)*exp(y(15))/exp(y(1))+exp(y(44))*(exp(y(21))-exp(y(41)));
T157 = exp(y(43))*(exp(y(1))*exp(y(44))*exp(y(35)))^params(7);
T169 = (y(42)+params(34))/(params(34)+y(9))-1;
T171 = params(9)/2*T169^2;
T182 = params(1)*exp(y(48))*params(9)*((params(34)+y(56))/(y(42)+params(34))-1);
T183 = ((params(34)+y(56))/(y(42)+params(34)))^2;
T184 = T182*T183;
T196 = exp(y(34))*params(7)*exp(y(15))/exp(y(35));
T200 = exp(y(1))*exp(y(44))*params(31)*exp(y(35))^params(5);
T238 = params(11)*exp(y(7))*exp(y(13))^((-params(12))*params(10));
T242 = T238*exp(y(46))^params(10);
T251 = (1-params(11)*exp(y(13))^(params(12)*(1-params(11)))*exp(y(46))^(params(11)-1))/(1-params(11));
T276 = exp(y(48))*params(1)*params(11)*exp(y(57))^params(10)*exp(y(46))^(params(12)*(-params(10)))*exp(y(54));
T283 = exp(y(48))*params(1)*params(11)*exp(y(57))^(params(10)-1);
T290 = T283*exp(y(46))^(params(12)*(1-params(10)))*exp(y(55));
T298 = exp(y(46))*exp(y(38))*params(10)/(params(10)-1)/exp(y(39));
T314 = exp(y(8))^params(15);
T318 = 1/params(1)*exp(y(46))^params(13);
T319 = exp(y(36))/(params(10)/(params(10)-1));
T324 = (T318*T319^params(14))^(1-params(15));
T352 = (-(exp(y(34))*params(7)*exp(y(15))/exp(y(1))/exp(y(2))));
T364 = (-(exp(y(17))^(1-params(7))*exp(y(43))*exp(y(1))*exp(y(44))*exp(y(35))*getPowerDeriv(exp(y(1))*exp(y(44))*exp(y(35)),params(7),1)));
T469 = getPowerDeriv(T318*T319^params(14),1-params(15),1);
T496 = params(9)/2*(-(y(42)+params(34)))/((params(34)+y(9))*(params(34)+y(9)))*2*T169;
T508 = params(9)/2*2*T169*1/(params(34)+y(9));
T552 = getPowerDeriv(T251,(-params(10))/(1-params(11)),1);
lhs =T18;
rhs =T26;
residual(1)= lhs-rhs;
lhs =exp(y(22));
rhs =params(35)*exp(y(19))^(params(35)-1)*(1-exp(y(17)))^(1-params(35));
residual(2)= lhs-rhs;
lhs =params(1)*exp(y(25))*exp(y(48));
rhs =1;
residual(3)= lhs-rhs;
lhs =exp(y(23));
rhs =exp(y(22))/exp(y(3));
residual(4)= lhs-rhs;
lhs =exp(y(29));
rhs =T68;
residual(5)= lhs-rhs;
lhs =exp(y(30));
rhs =1-params(6)+params(6)*params(1)*exp(y(48))*exp(y(52))*exp(y(51));
residual(6)= lhs-rhs;
lhs =exp(y(31));
rhs =exp(y(30))/(params(29)-exp(y(29)));
residual(7)= lhs-rhs;
lhs =exp(y(32));
rhs =exp(y(4))+(exp(y(24))-exp(y(4)))*exp(y(6));
residual(8)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(32))*exp(y(31))/exp(y(6));
residual(9)= lhs-rhs;
lhs =exp(y(21))*exp(y(16));
rhs =exp(y(31))*exp(y(26));
residual(10)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(27))+exp(y(28));
residual(11)= lhs-rhs;
lhs =exp(y(27));
rhs =params(6)*exp(y(32))*exp(y(5))*exp((-x(it_, 4)));
residual(12)= lhs-rhs;
lhs =exp(y(28));
rhs =exp(y(21))*params(28)*exp(y(44))*exp(y(1));
residual(13)= lhs-rhs;
lhs =exp(y(24));
rhs =T145/exp(y(2));
residual(14)= lhs-rhs;
lhs =exp(y(15));
rhs =T157*exp(y(17))^(1-params(7));
residual(15)= lhs-rhs;
lhs =exp(y(21));
rhs =1+T171+(y(42)+params(34))*params(9)*T169/(params(34)+y(9))-T184;
residual(16)= lhs-rhs;
lhs =exp(y(41));
rhs =params(32)+params(31)/(1+params(5))*exp(y(35))^(1+params(5));
residual(17)= lhs-rhs;
lhs =T196;
rhs =T200;
residual(18)= lhs-rhs;
lhs =y(42);
rhs =exp(y(18))-exp(y(1))*exp(y(44))*exp(y(41));
residual(19)= lhs-rhs;
lhs =exp(y(16));
rhs =y(42)+exp(y(44))*exp(y(1));
residual(20)= lhs-rhs;
lhs =exp(y(20));
rhs =params(33)*exp(y(45));
residual(21)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(18))+exp(y(19))+exp(y(20))+(y(42)+params(34))*T171;
residual(22)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(14))*exp(y(37));
residual(23)= lhs-rhs;
lhs =exp(y(37));
rhs =T242+(1-params(11))*T251^((-params(10))/(1-params(11)));
residual(24)= lhs-rhs;
lhs =exp(y(36));
rhs =1/exp(y(34));
residual(25)= lhs-rhs;
lhs =exp(y(38));
rhs =exp(y(34))*exp(y(14))+T276;
residual(26)= lhs-rhs;
lhs =exp(y(39));
rhs =exp(y(14))+T290;
residual(27)= lhs-rhs;
lhs =exp(y(47));
rhs =T298;
residual(28)= lhs-rhs;
lhs =exp(y(46))^(1-params(10));
rhs =params(11)*exp(y(13))^(params(12)*(1-params(10)))+(1-params(11))*exp(y(47))^(1-params(10));
residual(29)= lhs-rhs;
lhs =exp(y(40));
rhs =exp(y(25))*exp(y(57));
residual(30)= lhs-rhs;
lhs =exp(y(40));
rhs =T314*T324*exp(x(it_, 5));
residual(31)= lhs-rhs;
lhs =y(43);
rhs =params(18)*y(10)-x(it_, 1);
residual(32)= lhs-rhs;
lhs =y(44);
rhs =params(16)*y(11)-x(it_, 2);
residual(33)= lhs-rhs;
lhs =y(45);
rhs =params(20)*y(12)-x(it_, 3);
residual(34)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(34, 62);

  %
  % Jacobian matrix
  %

  g1(1,14)=T18;
  g1(1,17)=(-(exp(y(34))*(1-params(7))*exp(y(14))*exp(y(17))))/(exp(y(17))*exp(y(17)))-(1-params(35))*(-(exp(y(19))*(-exp(y(17)))))/((1-exp(y(17)))*(1-exp(y(17))))/params(35);
  g1(1,19)=(-T26);
  g1(1,34)=T18;
  g1(2,17)=(-(params(35)*exp(y(19))^(params(35)-1)*(-exp(y(17)))*getPowerDeriv(1-exp(y(17)),1-params(35),1)));
  g1(2,19)=(-((1-exp(y(17)))^(1-params(35))*params(35)*exp(y(19))*getPowerDeriv(exp(y(19)),params(35)-1,1)));
  g1(2,22)=exp(y(22));
  g1(3,48)=params(1)*exp(y(25))*exp(y(48));
  g1(3,25)=params(1)*exp(y(25))*exp(y(48));
  g1(4,3)=(-((-(exp(y(22))*exp(y(3))))/(exp(y(3))*exp(y(3)))));
  g1(4,22)=(-(exp(y(22))/exp(y(3))));
  g1(4,23)=exp(y(23));
  g1(5,48)=(-T68);
  g1(5,49)=(-(exp(y(48))*params(1)*(1-params(6))*exp(y(49))));
  g1(5,25)=(-(exp(y(48))*params(1)*(1-params(6))*(-exp(y(25)))));
  g1(5,29)=exp(y(29));
  g1(5,50)=(-(params(6)*params(1)*exp(y(48))*exp(y(53))*exp(y(50))));
  g1(5,53)=(-(params(6)*params(1)*exp(y(48))*exp(y(53))*exp(y(50))));
  g1(6,48)=(-(params(6)*params(1)*exp(y(48))*exp(y(52))*exp(y(51))));
  g1(6,30)=exp(y(30));
  g1(6,51)=(-(params(6)*params(1)*exp(y(48))*exp(y(52))*exp(y(51))));
  g1(6,52)=(-(params(6)*params(1)*exp(y(48))*exp(y(52))*exp(y(51))));
  g1(7,29)=(-((-(exp(y(30))*(-exp(y(29)))))/((params(29)-exp(y(29)))*(params(29)-exp(y(29))))));
  g1(7,30)=(-(exp(y(30))/(params(29)-exp(y(29)))));
  g1(7,31)=exp(y(31));
  g1(8,24)=(-(exp(y(24))*exp(y(6))));
  g1(8,4)=(-(exp(y(4))+exp(y(6))*(-exp(y(4)))));
  g1(8,6)=(-((exp(y(24))-exp(y(4)))*exp(y(6))));
  g1(8,32)=exp(y(32));
  g1(9,6)=(-(exp(y(32))*(-(exp(y(31))*exp(y(6))))/(exp(y(6))*exp(y(6)))));
  g1(9,31)=(-(exp(y(32))*exp(y(31))/exp(y(6))));
  g1(9,32)=(-(exp(y(32))*exp(y(31))/exp(y(6))));
  g1(9,33)=exp(y(33));
  g1(10,16)=exp(y(21))*exp(y(16));
  g1(10,21)=exp(y(21))*exp(y(16));
  g1(10,26)=(-(exp(y(31))*exp(y(26))));
  g1(10,31)=(-(exp(y(31))*exp(y(26))));
  g1(11,26)=exp(y(26));
  g1(11,27)=(-exp(y(27)));
  g1(11,28)=(-exp(y(28)));
  g1(12,5)=(-(params(6)*exp(y(32))*exp(y(5))*exp((-x(it_, 4)))));
  g1(12,27)=exp(y(27));
  g1(12,32)=(-(params(6)*exp(y(32))*exp(y(5))*exp((-x(it_, 4)))));
  g1(12,61)=(-(params(6)*exp(y(32))*exp(y(5))*(-exp((-x(it_, 4))))));
  g1(13,1)=(-(exp(y(21))*params(28)*exp(y(44))*exp(y(1))));
  g1(13,21)=(-(exp(y(21))*params(28)*exp(y(44))*exp(y(1))));
  g1(13,28)=exp(y(28));
  g1(13,44)=(-(exp(y(21))*params(28)*exp(y(44))*exp(y(1))));
  g1(14,15)=T352;
  g1(14,1)=(-((-(exp(y(1))*exp(y(34))*params(7)*exp(y(15))))/(exp(y(1))*exp(y(1)))/exp(y(2))));
  g1(14,2)=(-((-(T145*exp(y(2))))/(exp(y(2))*exp(y(2)))));
  g1(14,21)=(-(exp(y(21))*exp(y(44))/exp(y(2))));
  g1(14,24)=exp(y(24));
  g1(14,34)=T352;
  g1(14,41)=(-(exp(y(44))*(-exp(y(41)))/exp(y(2))));
  g1(14,44)=(-(exp(y(44))*(exp(y(21))-exp(y(41)))/exp(y(2))));
  g1(15,15)=exp(y(15));
  g1(15,1)=T364;
  g1(15,17)=(-(T157*exp(y(17))*getPowerDeriv(exp(y(17)),1-params(7),1)));
  g1(15,35)=T364;
  g1(15,43)=(-(T157*exp(y(17))^(1-params(7))));
  g1(15,44)=T364;
  g1(16,21)=exp(y(21));
  g1(16,48)=T184;
  g1(16,9)=(-(T496+((params(34)+y(9))*(y(42)+params(34))*params(9)*(-(y(42)+params(34)))/((params(34)+y(9))*(params(34)+y(9)))-(y(42)+params(34))*params(9)*T169)/((params(34)+y(9))*(params(34)+y(9)))));
  g1(16,42)=(-(T508+(params(9)*T169+(y(42)+params(34))*params(9)*1/(params(34)+y(9)))/(params(34)+y(9))-(T183*params(1)*exp(y(48))*params(9)*(-(params(34)+y(56)))/((y(42)+params(34))*(y(42)+params(34)))+T182*(-(params(34)+y(56)))/((y(42)+params(34))*(y(42)+params(34)))*2*(params(34)+y(56))/(y(42)+params(34)))));
  g1(16,56)=T183*params(1)*exp(y(48))*params(9)*1/(y(42)+params(34))+T182*2*(params(34)+y(56))/(y(42)+params(34))*1/(y(42)+params(34));
  g1(17,35)=(-(params(31)/(1+params(5))*exp(y(35))*getPowerDeriv(exp(y(35)),1+params(5),1)));
  g1(17,41)=exp(y(41));
  g1(18,15)=T196;
  g1(18,1)=(-T200);
  g1(18,34)=T196;
  g1(18,35)=(-(exp(y(34))*params(7)*exp(y(15))*exp(y(35))))/(exp(y(35))*exp(y(35)))-exp(y(1))*exp(y(44))*params(31)*exp(y(35))*getPowerDeriv(exp(y(35)),params(5),1);
  g1(18,44)=(-T200);
  g1(19,1)=exp(y(1))*exp(y(44))*exp(y(41));
  g1(19,18)=(-exp(y(18)));
  g1(19,41)=exp(y(1))*exp(y(44))*exp(y(41));
  g1(19,42)=1;
  g1(19,44)=exp(y(1))*exp(y(44))*exp(y(41));
  g1(20,1)=(-(exp(y(44))*exp(y(1))));
  g1(20,16)=exp(y(16));
  g1(20,42)=(-1);
  g1(20,44)=(-(exp(y(44))*exp(y(1))));
  g1(21,20)=exp(y(20));
  g1(21,45)=(-(params(33)*exp(y(45))));
  g1(22,14)=exp(y(14));
  g1(22,18)=(-exp(y(18)));
  g1(22,19)=(-exp(y(19)));
  g1(22,20)=(-exp(y(20)));
  g1(22,9)=(-((y(42)+params(34))*T496));
  g1(22,42)=(-(T171+(y(42)+params(34))*T508));
  g1(23,14)=(-(exp(y(14))*exp(y(37))));
  g1(23,15)=exp(y(15));
  g1(23,37)=(-(exp(y(14))*exp(y(37))));
  g1(24,7)=(-T242);
  g1(24,37)=exp(y(37));
  g1(24,13)=(-(exp(y(46))^params(10)*params(11)*exp(y(7))*exp(y(13))*getPowerDeriv(exp(y(13)),(-params(12))*params(10),1)+(1-params(11))*(-(exp(y(46))^(params(11)-1)*params(11)*exp(y(13))*getPowerDeriv(exp(y(13)),params(12)*(1-params(11)),1)))/(1-params(11))*T552));
  g1(24,46)=(-(T238*exp(y(46))*getPowerDeriv(exp(y(46)),params(10),1)+(1-params(11))*T552*(-(params(11)*exp(y(13))^(params(12)*(1-params(11)))*exp(y(46))*getPowerDeriv(exp(y(46)),params(11)-1,1)))/(1-params(11))));
  g1(25,34)=(-((-exp(y(34)))/(exp(y(34))*exp(y(34)))));
  g1(25,36)=exp(y(36));
  g1(26,14)=(-(exp(y(34))*exp(y(14))));
  g1(26,48)=(-T276);
  g1(26,34)=(-(exp(y(34))*exp(y(14))));
  g1(26,38)=exp(y(38));
  g1(26,54)=(-T276);
  g1(26,46)=(-(exp(y(54))*exp(y(48))*params(1)*params(11)*exp(y(57))^params(10)*exp(y(46))*getPowerDeriv(exp(y(46)),params(12)*(-params(10)),1)));
  g1(26,57)=(-(exp(y(54))*exp(y(46))^(params(12)*(-params(10)))*exp(y(48))*params(1)*params(11)*exp(y(57))*getPowerDeriv(exp(y(57)),params(10),1)));
  g1(27,14)=(-exp(y(14)));
  g1(27,48)=(-T290);
  g1(27,39)=exp(y(39));
  g1(27,55)=(-T290);
  g1(27,46)=(-(exp(y(55))*T283*exp(y(46))*getPowerDeriv(exp(y(46)),params(12)*(1-params(10)),1)));
  g1(27,57)=(-(exp(y(55))*exp(y(46))^(params(12)*(1-params(10)))*exp(y(48))*params(1)*params(11)*exp(y(57))*getPowerDeriv(exp(y(57)),params(10)-1,1)));
  g1(28,38)=(-T298);
  g1(28,39)=(-(exp(y(46))*(-(exp(y(39))*exp(y(38))*params(10)/(params(10)-1)))/(exp(y(39))*exp(y(39)))));
  g1(28,46)=(-T298);
  g1(28,47)=exp(y(47));
  g1(29,13)=(-(params(11)*exp(y(13))*getPowerDeriv(exp(y(13)),params(12)*(1-params(10)),1)));
  g1(29,46)=exp(y(46))*getPowerDeriv(exp(y(46)),1-params(10),1);
  g1(29,47)=(-((1-params(11))*exp(y(47))*getPowerDeriv(exp(y(47)),1-params(10),1)));
  g1(30,25)=(-(exp(y(25))*exp(y(57))));
  g1(30,40)=exp(y(40));
  g1(30,57)=(-(exp(y(25))*exp(y(57))));
  g1(31,36)=(-(exp(x(it_, 5))*T314*T318*T319*getPowerDeriv(T319,params(14),1)*T469));
  g1(31,8)=(-(exp(x(it_, 5))*T324*exp(y(8))*getPowerDeriv(exp(y(8)),params(15),1)));
  g1(31,40)=exp(y(40));
  g1(31,46)=(-(exp(x(it_, 5))*T314*T469*T319^params(14)*1/params(1)*exp(y(46))*getPowerDeriv(exp(y(46)),params(13),1)));
  g1(31,62)=(-(T314*T324*exp(x(it_, 5))));
  g1(32,10)=(-params(18));
  g1(32,43)=1;
  g1(32,58)=1;
  g1(33,11)=(-params(16));
  g1(33,44)=1;
  g1(33,59)=1;
  g1(34,12)=(-params(20));
  g1(34,45)=1;
  g1(34,60)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],34,3844);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],34,238328);
end
end
