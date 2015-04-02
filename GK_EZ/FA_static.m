function [residual, g1, g2] = FA_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                    columns: equations in order of declaration
%                                                    rows: variables in declaration order
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: equations in order of declaration
%                                                       rows: variables in declaration order
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 36, 1);

%
% Model equations
%

T26 = y(34)^(1-params(36));
T37 = (1-params(1))*y(35)^((1-params(36))/params(37))+params(1)*y(36)^(1/params(37));
T43 = T26/y(36);
T54 = exp(y(6))*(1-params(34))/params(34)/(1-exp(y(4)));
T80 = exp(y(9))*params(1)*(1-params(5))*(exp(y(10))-exp(y(11)))+exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19));
T138 = exp(y(20))*params(6)*exp(y(2))/exp(y(3))+exp(y(30))*(exp(y(8))-exp(y(27)));
T148 = exp(y(29))*(exp(y(3))*exp(y(30))*exp(y(21)))^params(6);
T168 = exp(y(20))*params(6)*exp(y(2))/exp(y(21));
T172 = exp(y(3))*exp(y(30))*params(30)*exp(y(21))^params(4);
T200 = exp(y(32));
T206 = exp(y(23))*params(10)*T200^((-params(11))*params(9));
T207 = T200^params(9);
T208 = T206*T207;
T217 = (1-params(10)*T200^(params(11)*(1-params(10)))*T200^(params(10)-1))/(1-params(10));
T237 = exp(y(24))*T207*exp(y(9))*params(1)*params(10)*T200^(params(11)*(-params(9)));
T247 = T200^(params(11)*(1-params(9)));
T249 = exp(y(25))*exp(y(9))*params(1)*params(10)*T200^(params(9)-1)*T247;
T257 = T200*exp(y(24))*params(9)/(params(9)-1)/exp(y(25));
T270 = exp(y(26))^params(14);
T274 = 1/params(1)*T200^params(12);
T275 = exp(y(22))/(params(9)/(params(9)-1));
T280 = (T274*T275^params(13))^(1-params(14));
T301 = (-(exp(y(20))*(1-params(6))*exp(y(1))/exp(y(4))));
T306 = (-(exp(y(20))*params(6)*exp(y(2))/exp(y(3))/exp(y(8))));
T318 = (-(exp(y(4))^(1-params(6))*exp(y(29))*exp(y(3))*exp(y(30))*exp(y(21))*getPowerDeriv(exp(y(3))*exp(y(30))*exp(y(21)),params(6),1)));
T407 = getPowerDeriv(T274*T275^params(13),1-params(14),1);
T499 = getPowerDeriv(T43,1-1/params(37),1);
T504 = getPowerDeriv(T37,params(37)/(1-params(36)),1);
lhs =y(35);
rhs =params(38)*exp(y(6))^params(34)*(1-exp(y(4)))^(1-params(34));
residual(1)= lhs-rhs;
lhs =y(36);
rhs =T26;
residual(2)= lhs-rhs;
lhs =y(34);
rhs =T37^(params(37)/(1-params(36)));
residual(3)= lhs-rhs;
lhs =exp(y(9));
rhs =T43^(1-1/params(37));
residual(4)= lhs-rhs;
lhs =params(1)*exp(y(9))*exp(y(11));
rhs =1;
residual(5)= lhs-rhs;
lhs =T54;
rhs =exp(y(20))*(1-params(6))*exp(y(1))/exp(y(4));
residual(6)= lhs-rhs;
lhs =exp(y(15));
rhs =T80;
residual(7)= lhs-rhs;
lhs =exp(y(16));
rhs =1-params(5)+exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18));
residual(8)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(16))/(params(28)-exp(y(15)));
residual(9)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(11))+(exp(y(10))-exp(y(11)))*exp(y(17));
residual(10)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(18));
residual(11)= lhs-rhs;
lhs =exp(y(8))*exp(y(3));
rhs =exp(y(17))*exp(y(12));
residual(12)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(13))+exp(y(14));
residual(13)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(12))*params(5)*exp(y(18))*exp((-x(4)));
residual(14)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(3))*exp(y(8))*params(27)*exp(y(30));
residual(15)= lhs-rhs;
lhs =exp(y(10));
rhs =T138/exp(y(8));
residual(16)= lhs-rhs;
lhs =exp(y(2));
rhs =T148*exp(y(4))^(1-params(6));
residual(17)= lhs-rhs;
lhs =exp(y(8));
rhs =1;
residual(18)= lhs-rhs;
lhs =exp(y(27));
rhs =params(31)+params(30)/(1+params(4))*exp(y(21))^(1+params(4));
residual(19)= lhs-rhs;
lhs =T168;
rhs =T172;
residual(20)= lhs-rhs;
lhs =y(28);
rhs =exp(y(5))-exp(y(3))*exp(y(30))*exp(y(27));
residual(21)= lhs-rhs;
lhs =exp(y(3));
rhs =y(28)+exp(y(3))*exp(y(30));
residual(22)= lhs-rhs;
lhs =exp(y(7));
rhs =params(32)*exp(y(31));
residual(23)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(5))+exp(y(6))+exp(y(7));
residual(24)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(1))*exp(y(23));
residual(25)= lhs-rhs;
lhs =exp(y(23));
rhs =T208+(1-params(10))*T217^((-params(9))/(1-params(10)));
residual(26)= lhs-rhs;
lhs =exp(y(22));
rhs =1/exp(y(20));
residual(27)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(20))*exp(y(1))+T237;
residual(28)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(1))+T249;
residual(29)= lhs-rhs;
lhs =exp(y(33));
rhs =T257;
residual(30)= lhs-rhs;
lhs =T200^(1-params(9));
rhs =params(10)*T247+(1-params(10))*exp(y(33))^(1-params(9));
residual(31)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(11))*T200;
residual(32)= lhs-rhs;
lhs =exp(y(26));
rhs =T270*T280*exp(x(5));
residual(33)= lhs-rhs;
lhs =y(29);
rhs =y(29)*params(17)-x(1);
residual(34)= lhs-rhs;
lhs =y(30);
rhs =y(30)*params(15)-x(2);
residual(35)= lhs-rhs;
lhs =y(31);
rhs =y(31)*params(19)-x(3);
residual(36)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(36, 36);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-(params(38)*exp(y(6))^params(34)*(-exp(y(4)))*getPowerDeriv(1-exp(y(4)),1-params(34),1)));
  g1(1,6)=(-((1-exp(y(4)))^(1-params(34))*params(38)*exp(y(6))*getPowerDeriv(exp(y(6)),params(34),1)));
  g1(1,35)=1;
  g1(2,34)=(-(getPowerDeriv(y(34),1-params(36),1)));
  g1(2,36)=1;
  g1(3,34)=1;
  g1(3,35)=(-((1-params(1))*getPowerDeriv(y(35),(1-params(36))/params(37),1)*T504));
  g1(3,36)=(-(T504*params(1)*getPowerDeriv(y(36),1/params(37),1)));
  g1(4,9)=exp(y(9));
  g1(4,34)=(-(getPowerDeriv(y(34),1-params(36),1)/y(36)*T499));
  g1(4,36)=(-(T499*(-T26)/(y(36)*y(36))));
  g1(5,9)=params(1)*exp(y(9))*exp(y(11));
  g1(5,11)=params(1)*exp(y(9))*exp(y(11));
  g1(6,1)=T301;
  g1(6,4)=(-(exp(y(6))*(1-params(34))/params(34)*(-exp(y(4)))))/((1-exp(y(4)))*(1-exp(y(4))))-(-(exp(y(4))*exp(y(20))*(1-params(6))*exp(y(1))))/(exp(y(4))*exp(y(4)));
  g1(6,6)=T54;
  g1(6,20)=T301;
  g1(7,9)=(-T80);
  g1(7,10)=(-(exp(y(9))*params(1)*(1-params(5))*exp(y(10))));
  g1(7,11)=(-(exp(y(9))*params(1)*(1-params(5))*(-exp(y(11)))));
  g1(7,15)=exp(y(15))-exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19));
  g1(7,19)=(-(exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19))));
  g1(8,9)=(-(exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18))));
  g1(8,16)=exp(y(16))-exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18));
  g1(8,18)=(-(exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18))));
  g1(9,15)=(-((-(exp(y(16))*(-exp(y(15)))))/((params(28)-exp(y(15)))*(params(28)-exp(y(15))))));
  g1(9,16)=(-(exp(y(16))/(params(28)-exp(y(15)))));
  g1(9,17)=exp(y(17));
  g1(10,10)=(-(exp(y(10))*exp(y(17))));
  g1(10,11)=(-(exp(y(11))+exp(y(17))*(-exp(y(11)))));
  g1(10,17)=(-((exp(y(10))-exp(y(11)))*exp(y(17))));
  g1(10,18)=exp(y(18));
  g1(11,18)=(-exp(y(18)));
  g1(11,19)=exp(y(19));
  g1(12,3)=exp(y(8))*exp(y(3));
  g1(12,8)=exp(y(8))*exp(y(3));
  g1(12,12)=(-(exp(y(17))*exp(y(12))));
  g1(12,17)=(-(exp(y(17))*exp(y(12))));
  g1(13,12)=exp(y(12));
  g1(13,13)=(-exp(y(13)));
  g1(13,14)=(-exp(y(14)));
  g1(14,12)=(-(exp(y(12))*params(5)*exp(y(18))*exp((-x(4)))));
  g1(14,13)=exp(y(13));
  g1(14,18)=(-(exp(y(12))*params(5)*exp(y(18))*exp((-x(4)))));
  g1(15,3)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(15,8)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(15,14)=exp(y(14));
  g1(15,30)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(16,2)=T306;
  g1(16,3)=(-((-(exp(y(3))*exp(y(20))*params(6)*exp(y(2))))/(exp(y(3))*exp(y(3)))/exp(y(8))));
  g1(16,8)=(-((exp(y(8))*exp(y(8))*exp(y(30))-exp(y(8))*T138)/(exp(y(8))*exp(y(8)))));
  g1(16,10)=exp(y(10));
  g1(16,20)=T306;
  g1(16,27)=(-(exp(y(30))*(-exp(y(27)))/exp(y(8))));
  g1(16,30)=(-(exp(y(30))*(exp(y(8))-exp(y(27)))/exp(y(8))));
  g1(17,2)=exp(y(2));
  g1(17,3)=T318;
  g1(17,4)=(-(T148*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(6),1)));
  g1(17,21)=T318;
  g1(17,29)=(-(T148*exp(y(4))^(1-params(6))));
  g1(17,30)=T318;
  g1(18,8)=exp(y(8));
  g1(19,21)=(-(params(30)/(1+params(4))*exp(y(21))*getPowerDeriv(exp(y(21)),1+params(4),1)));
  g1(19,27)=exp(y(27));
  g1(20,2)=T168;
  g1(20,3)=(-T172);
  g1(20,20)=T168;
  g1(20,21)=(-(exp(y(20))*params(6)*exp(y(2))*exp(y(21))))/(exp(y(21))*exp(y(21)))-exp(y(3))*exp(y(30))*params(30)*exp(y(21))*getPowerDeriv(exp(y(21)),params(4),1);
  g1(20,30)=(-T172);
  g1(21,3)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(21,5)=(-exp(y(5)));
  g1(21,27)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(21,28)=1;
  g1(21,30)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(22,3)=exp(y(3))-exp(y(3))*exp(y(30));
  g1(22,28)=(-1);
  g1(22,30)=(-(exp(y(3))*exp(y(30))));
  g1(23,7)=exp(y(7));
  g1(23,31)=(-(params(32)*exp(y(31))));
  g1(24,1)=exp(y(1));
  g1(24,5)=(-exp(y(5)));
  g1(24,6)=(-exp(y(6)));
  g1(24,7)=(-exp(y(7)));
  g1(25,1)=(-(exp(y(1))*exp(y(23))));
  g1(25,2)=exp(y(2));
  g1(25,23)=(-(exp(y(1))*exp(y(23))));
  g1(26,23)=exp(y(23))-T208;
  g1(26,32)=(-(T207*exp(y(23))*params(10)*T200*getPowerDeriv(T200,(-params(11))*params(9),1)+T206*T200*getPowerDeriv(T200,params(9),1)+(1-params(10))*(-(T200^(params(10)-1)*params(10)*T200*getPowerDeriv(T200,params(11)*(1-params(10)),1)+params(10)*T200^(params(11)*(1-params(10)))*T200*getPowerDeriv(T200,params(10)-1,1)))/(1-params(10))*getPowerDeriv(T217,(-params(9))/(1-params(10)),1)));
  g1(27,20)=(-((-exp(y(20)))/(exp(y(20))*exp(y(20)))));
  g1(27,22)=exp(y(22));
  g1(28,1)=(-(exp(y(20))*exp(y(1))));
  g1(28,9)=(-T237);
  g1(28,20)=(-(exp(y(20))*exp(y(1))));
  g1(28,24)=exp(y(24))-T237;
  g1(28,32)=(-(exp(y(24))*(T200^(params(11)*(-params(9)))*exp(y(9))*params(1)*params(10)*T200*getPowerDeriv(T200,params(9),1)+T207*exp(y(9))*params(1)*params(10)*T200*getPowerDeriv(T200,params(11)*(-params(9)),1))));
  g1(29,1)=(-exp(y(1)));
  g1(29,9)=(-T249);
  g1(29,25)=exp(y(25))-T249;
  g1(29,32)=(-(exp(y(25))*(T247*exp(y(9))*params(1)*params(10)*T200*getPowerDeriv(T200,params(9)-1,1)+exp(y(9))*params(1)*params(10)*T200^(params(9)-1)*T200*getPowerDeriv(T200,params(11)*(1-params(9)),1))));
  g1(30,24)=(-T257);
  g1(30,25)=(-(T200*(-(exp(y(25))*exp(y(24))*params(9)/(params(9)-1)))/(exp(y(25))*exp(y(25)))));
  g1(30,32)=(-T257);
  g1(30,33)=exp(y(33));
  g1(31,32)=T200*getPowerDeriv(T200,1-params(9),1)-params(10)*T200*getPowerDeriv(T200,params(11)*(1-params(9)),1);
  g1(31,33)=(-((1-params(10))*exp(y(33))*getPowerDeriv(exp(y(33)),1-params(9),1)));
  g1(32,11)=(-(exp(y(11))*T200));
  g1(32,26)=exp(y(26));
  g1(32,32)=(-(exp(y(11))*T200));
  g1(33,22)=(-(exp(x(5))*T270*T274*T275*getPowerDeriv(T275,params(13),1)*T407));
  g1(33,26)=exp(y(26))-exp(x(5))*T280*exp(y(26))*getPowerDeriv(exp(y(26)),params(14),1);
  g1(33,32)=(-(exp(x(5))*T270*T407*T275^params(13)*1/params(1)*T200*getPowerDeriv(T200,params(12),1)));
  g1(34,29)=1-params(17);
  g1(35,30)=1-params(15);
  g1(36,31)=1-params(19);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],36,1296);
end
end
