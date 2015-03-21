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

residual = zeros( 34, 1);

%
% Model equations
%

T20 = exp(y(6))^params(34)*(1-exp(y(4)))^(1-params(34));
T31 = (1-params(1))*T20^((1-params(36))/params(37))+params(1)*(y(34)^(1-params(36)))^(1/params(37));
T37 = exp(y(6))*(1-params(34))/params(34)/(1-exp(y(4)));
T71 = exp(y(9))*params(1)*(1-params(5))*(exp(y(10))-exp(y(11)))+exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19));
T129 = exp(y(20))*params(6)*exp(y(2))/exp(y(3))+exp(y(30))*(exp(y(8))-exp(y(27)));
T139 = exp(y(29))*(exp(y(3))*exp(y(30))*exp(y(21)))^params(6);
T159 = exp(y(20))*params(6)*exp(y(2))/exp(y(21));
T163 = exp(y(3))*exp(y(30))*params(30)*exp(y(21))^params(4);
T191 = exp(y(32));
T197 = exp(y(23))*params(10)*T191^((-params(11))*params(9));
T198 = T191^params(9);
T199 = T197*T198;
T208 = (1-params(10)*T191^(params(11)*(1-params(10)))*T191^(params(10)-1))/(1-params(10));
T228 = exp(y(24))*T198*exp(y(9))*params(1)*params(10)*T191^(params(11)*(-params(9)));
T238 = T191^(params(11)*(1-params(9)));
T240 = exp(y(25))*exp(y(9))*params(1)*params(10)*T191^(params(9)-1)*T238;
T248 = T191*exp(y(24))*params(9)/(params(9)-1)/exp(y(25));
T261 = exp(y(26))^params(14);
T265 = 1/params(1)*T191^params(12);
T266 = exp(y(22))/(params(9)/(params(9)-1));
T271 = (T265*T266^params(13))^(1-params(14));
T292 = (-(exp(y(20))*(1-params(6))*exp(y(1))/exp(y(4))));
T297 = (-(exp(y(20))*params(6)*exp(y(2))/exp(y(3))/exp(y(8))));
T309 = (-(exp(y(4))^(1-params(6))*exp(y(29))*exp(y(3))*exp(y(30))*exp(y(21))*getPowerDeriv(exp(y(3))*exp(y(30))*exp(y(21)),params(6),1)));
T317 = getPowerDeriv(T20,(1-params(36))/params(37),1);
T320 = getPowerDeriv(T31,params(37)/(1-params(36)),1);
T405 = getPowerDeriv(T265*T266^params(13),1-params(14),1);
lhs =y(34);
rhs =T31^(params(37)/(1-params(36)));
residual(1)= lhs-rhs;
lhs =T37;
rhs =exp(y(20))*(1-params(6))*exp(y(1))/exp(y(4));
residual(2)= lhs-rhs;
lhs =exp(y(9));
rhs =1;
residual(3)= lhs-rhs;
lhs =params(1)*exp(y(9))*exp(y(11));
rhs =1;
residual(4)= lhs-rhs;
lhs =exp(y(15));
rhs =T71;
residual(5)= lhs-rhs;
lhs =exp(y(16));
rhs =1-params(5)+exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18));
residual(6)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(16))/(params(28)-exp(y(15)));
residual(7)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(11))+(exp(y(10))-exp(y(11)))*exp(y(17));
residual(8)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(18));
residual(9)= lhs-rhs;
lhs =exp(y(8))*exp(y(3));
rhs =exp(y(17))*exp(y(12));
residual(10)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(13))+exp(y(14));
residual(11)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(12))*params(5)*exp(y(18))*exp((-x(4)));
residual(12)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(3))*exp(y(8))*params(27)*exp(y(30));
residual(13)= lhs-rhs;
lhs =exp(y(10));
rhs =T129/exp(y(8));
residual(14)= lhs-rhs;
lhs =exp(y(2));
rhs =T139*exp(y(4))^(1-params(6));
residual(15)= lhs-rhs;
lhs =exp(y(8));
rhs =1;
residual(16)= lhs-rhs;
lhs =exp(y(27));
rhs =params(31)+params(30)/(1+params(4))*exp(y(21))^(1+params(4));
residual(17)= lhs-rhs;
lhs =T159;
rhs =T163;
residual(18)= lhs-rhs;
lhs =y(28);
rhs =exp(y(5))-exp(y(3))*exp(y(30))*exp(y(27));
residual(19)= lhs-rhs;
lhs =exp(y(3));
rhs =y(28)+exp(y(3))*exp(y(30));
residual(20)= lhs-rhs;
lhs =exp(y(7));
rhs =params(32)*exp(y(31));
residual(21)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(5))+exp(y(6))+exp(y(7));
residual(22)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(1))*exp(y(23));
residual(23)= lhs-rhs;
lhs =exp(y(23));
rhs =T199+(1-params(10))*T208^((-params(9))/(1-params(10)));
residual(24)= lhs-rhs;
lhs =exp(y(22));
rhs =1/exp(y(20));
residual(25)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(20))*exp(y(1))+T228;
residual(26)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(1))+T240;
residual(27)= lhs-rhs;
lhs =exp(y(33));
rhs =T248;
residual(28)= lhs-rhs;
lhs =T191^(1-params(9));
rhs =params(10)*T238+(1-params(10))*exp(y(33))^(1-params(9));
residual(29)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(11))*T191;
residual(30)= lhs-rhs;
lhs =exp(y(26));
rhs =T261*T271*exp(x(5));
residual(31)= lhs-rhs;
lhs =y(29);
rhs =y(29)*params(17)-x(1);
residual(32)= lhs-rhs;
lhs =y(30);
rhs =y(30)*params(15)-x(2);
residual(33)= lhs-rhs;
lhs =y(31);
rhs =y(31)*params(19)-x(3);
residual(34)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(34, 34);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-((1-params(1))*exp(y(6))^params(34)*(-exp(y(4)))*getPowerDeriv(1-exp(y(4)),1-params(34),1)*T317*T320));
  g1(1,6)=(-(T320*(1-params(1))*T317*(1-exp(y(4)))^(1-params(34))*exp(y(6))*getPowerDeriv(exp(y(6)),params(34),1)));
  g1(1,34)=1-T320*params(1)*getPowerDeriv(y(34),1-params(36),1)*getPowerDeriv(y(34)^(1-params(36)),1/params(37),1);
  g1(2,1)=T292;
  g1(2,4)=(-(exp(y(6))*(1-params(34))/params(34)*(-exp(y(4)))))/((1-exp(y(4)))*(1-exp(y(4))))-(-(exp(y(4))*exp(y(20))*(1-params(6))*exp(y(1))))/(exp(y(4))*exp(y(4)));
  g1(2,6)=T37;
  g1(2,20)=T292;
  g1(3,9)=exp(y(9));
  g1(4,9)=params(1)*exp(y(9))*exp(y(11));
  g1(4,11)=params(1)*exp(y(9))*exp(y(11));
  g1(5,9)=(-T71);
  g1(5,10)=(-(exp(y(9))*params(1)*(1-params(5))*exp(y(10))));
  g1(5,11)=(-(exp(y(9))*params(1)*(1-params(5))*(-exp(y(11)))));
  g1(5,15)=exp(y(15))-exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19));
  g1(5,19)=(-(exp(y(15))*params(1)*exp(y(9))*params(5)*exp(y(19))));
  g1(6,9)=(-(exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18))));
  g1(6,16)=exp(y(16))-exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18));
  g1(6,18)=(-(exp(y(16))*params(1)*exp(y(9))*params(5)*exp(y(18))));
  g1(7,15)=(-((-(exp(y(16))*(-exp(y(15)))))/((params(28)-exp(y(15)))*(params(28)-exp(y(15))))));
  g1(7,16)=(-(exp(y(16))/(params(28)-exp(y(15)))));
  g1(7,17)=exp(y(17));
  g1(8,10)=(-(exp(y(10))*exp(y(17))));
  g1(8,11)=(-(exp(y(11))+exp(y(17))*(-exp(y(11)))));
  g1(8,17)=(-((exp(y(10))-exp(y(11)))*exp(y(17))));
  g1(8,18)=exp(y(18));
  g1(9,18)=(-exp(y(18)));
  g1(9,19)=exp(y(19));
  g1(10,3)=exp(y(8))*exp(y(3));
  g1(10,8)=exp(y(8))*exp(y(3));
  g1(10,12)=(-(exp(y(17))*exp(y(12))));
  g1(10,17)=(-(exp(y(17))*exp(y(12))));
  g1(11,12)=exp(y(12));
  g1(11,13)=(-exp(y(13)));
  g1(11,14)=(-exp(y(14)));
  g1(12,12)=(-(exp(y(12))*params(5)*exp(y(18))*exp((-x(4)))));
  g1(12,13)=exp(y(13));
  g1(12,18)=(-(exp(y(12))*params(5)*exp(y(18))*exp((-x(4)))));
  g1(13,3)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(13,8)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(13,14)=exp(y(14));
  g1(13,30)=(-(exp(y(3))*exp(y(8))*params(27)*exp(y(30))));
  g1(14,2)=T297;
  g1(14,3)=(-((-(exp(y(3))*exp(y(20))*params(6)*exp(y(2))))/(exp(y(3))*exp(y(3)))/exp(y(8))));
  g1(14,8)=(-((exp(y(8))*exp(y(8))*exp(y(30))-exp(y(8))*T129)/(exp(y(8))*exp(y(8)))));
  g1(14,10)=exp(y(10));
  g1(14,20)=T297;
  g1(14,27)=(-(exp(y(30))*(-exp(y(27)))/exp(y(8))));
  g1(14,30)=(-(exp(y(30))*(exp(y(8))-exp(y(27)))/exp(y(8))));
  g1(15,2)=exp(y(2));
  g1(15,3)=T309;
  g1(15,4)=(-(T139*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(6),1)));
  g1(15,21)=T309;
  g1(15,29)=(-(T139*exp(y(4))^(1-params(6))));
  g1(15,30)=T309;
  g1(16,8)=exp(y(8));
  g1(17,21)=(-(params(30)/(1+params(4))*exp(y(21))*getPowerDeriv(exp(y(21)),1+params(4),1)));
  g1(17,27)=exp(y(27));
  g1(18,2)=T159;
  g1(18,3)=(-T163);
  g1(18,20)=T159;
  g1(18,21)=(-(exp(y(20))*params(6)*exp(y(2))*exp(y(21))))/(exp(y(21))*exp(y(21)))-exp(y(3))*exp(y(30))*params(30)*exp(y(21))*getPowerDeriv(exp(y(21)),params(4),1);
  g1(18,30)=(-T163);
  g1(19,3)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(19,5)=(-exp(y(5)));
  g1(19,27)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(19,28)=1;
  g1(19,30)=exp(y(3))*exp(y(30))*exp(y(27));
  g1(20,3)=exp(y(3))-exp(y(3))*exp(y(30));
  g1(20,28)=(-1);
  g1(20,30)=(-(exp(y(3))*exp(y(30))));
  g1(21,7)=exp(y(7));
  g1(21,31)=(-(params(32)*exp(y(31))));
  g1(22,1)=exp(y(1));
  g1(22,5)=(-exp(y(5)));
  g1(22,6)=(-exp(y(6)));
  g1(22,7)=(-exp(y(7)));
  g1(23,1)=(-(exp(y(1))*exp(y(23))));
  g1(23,2)=exp(y(2));
  g1(23,23)=(-(exp(y(1))*exp(y(23))));
  g1(24,23)=exp(y(23))-T199;
  g1(24,32)=(-(T198*exp(y(23))*params(10)*T191*getPowerDeriv(T191,(-params(11))*params(9),1)+T197*T191*getPowerDeriv(T191,params(9),1)+(1-params(10))*(-(T191^(params(10)-1)*params(10)*T191*getPowerDeriv(T191,params(11)*(1-params(10)),1)+params(10)*T191^(params(11)*(1-params(10)))*T191*getPowerDeriv(T191,params(10)-1,1)))/(1-params(10))*getPowerDeriv(T208,(-params(9))/(1-params(10)),1)));
  g1(25,20)=(-((-exp(y(20)))/(exp(y(20))*exp(y(20)))));
  g1(25,22)=exp(y(22));
  g1(26,1)=(-(exp(y(20))*exp(y(1))));
  g1(26,9)=(-T228);
  g1(26,20)=(-(exp(y(20))*exp(y(1))));
  g1(26,24)=exp(y(24))-T228;
  g1(26,32)=(-(exp(y(24))*(T191^(params(11)*(-params(9)))*exp(y(9))*params(1)*params(10)*T191*getPowerDeriv(T191,params(9),1)+T198*exp(y(9))*params(1)*params(10)*T191*getPowerDeriv(T191,params(11)*(-params(9)),1))));
  g1(27,1)=(-exp(y(1)));
  g1(27,9)=(-T240);
  g1(27,25)=exp(y(25))-T240;
  g1(27,32)=(-(exp(y(25))*(T238*exp(y(9))*params(1)*params(10)*T191*getPowerDeriv(T191,params(9)-1,1)+exp(y(9))*params(1)*params(10)*T191^(params(9)-1)*T191*getPowerDeriv(T191,params(11)*(1-params(9)),1))));
  g1(28,24)=(-T248);
  g1(28,25)=(-(T191*(-(exp(y(25))*exp(y(24))*params(9)/(params(9)-1)))/(exp(y(25))*exp(y(25)))));
  g1(28,32)=(-T248);
  g1(28,33)=exp(y(33));
  g1(29,32)=T191*getPowerDeriv(T191,1-params(9),1)-params(10)*T191*getPowerDeriv(T191,params(11)*(1-params(9)),1);
  g1(29,33)=(-((1-params(10))*exp(y(33))*getPowerDeriv(exp(y(33)),1-params(9),1)));
  g1(30,11)=(-(exp(y(11))*T191));
  g1(30,26)=exp(y(26));
  g1(30,32)=(-(exp(y(11))*T191));
  g1(31,22)=(-(exp(x(5))*T261*T265*T266*getPowerDeriv(T266,params(13),1)*T405));
  g1(31,26)=exp(y(26))-exp(x(5))*T271*exp(y(26))*getPowerDeriv(exp(y(26)),params(14),1);
  g1(31,32)=(-(exp(x(5))*T261*T405*T266^params(13)*1/params(1)*T191*getPowerDeriv(T191,params(12),1)));
  g1(32,29)=1-params(17);
  g1(33,30)=1-params(15);
  g1(34,31)=1-params(19);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],34,1156);
end
end
