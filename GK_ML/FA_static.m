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

T18 = exp(y(21))*(1-params(7))*exp(y(1))/exp(y(4));
T26 = exp(y(6))/(1-exp(y(4)))*(1-params(35))/params(35);
T61 = exp(y(10))*params(1)*(1-params(6))*(exp(y(11))-exp(y(12)))+exp(y(16))*params(6)*params(1)*exp(y(10))*exp(y(20));
T119 = exp(y(21))*params(7)*exp(y(2))/exp(y(3))+exp(y(31))*(exp(y(8))-exp(y(28)));
T129 = exp(y(30))*(exp(y(3))*exp(y(31))*exp(y(22)))^params(7);
T149 = exp(y(21))*params(7)*exp(y(2))/exp(y(22));
T153 = exp(y(3))*exp(y(31))*params(31)*exp(y(22))^params(5);
T181 = exp(y(33));
T187 = exp(y(24))*params(11)*T181^((-params(12))*params(10));
T188 = T181^params(10);
T189 = T187*T188;
T198 = (1-params(11)*T181^(params(12)*(1-params(11)))*T181^(params(11)-1))/(1-params(11));
T218 = exp(y(25))*T188*exp(y(10))*params(1)*params(11)*T181^(params(12)*(-params(10)));
T228 = T181^(params(12)*(1-params(10)));
T230 = exp(y(26))*exp(y(10))*params(1)*params(11)*T181^(params(10)-1)*T228;
T238 = T181*exp(y(25))*params(10)/(params(10)-1)/exp(y(26));
T251 = exp(y(27))^params(15);
T255 = 1/params(1)*T181^params(13);
T256 = exp(y(23))/(params(10)/(params(10)-1));
T261 = (T255*T256^params(14))^(1-params(15));
T286 = (-(exp(y(21))*params(7)*exp(y(2))/exp(y(3))/exp(y(8))));
T298 = (-(exp(y(4))^(1-params(7))*exp(y(30))*exp(y(3))*exp(y(31))*exp(y(22))*getPowerDeriv(exp(y(3))*exp(y(31))*exp(y(22)),params(7),1)));
T390 = getPowerDeriv(T255*T256^params(14),1-params(15),1);
lhs =T18;
rhs =T26;
residual(1)= lhs-rhs;
lhs =exp(y(9));
rhs =params(35)*exp(y(6))^(params(35)-1)*(1-exp(y(4)))^(1-params(35));
residual(2)= lhs-rhs;
lhs =params(1)*exp(y(12))*exp(y(10));
rhs =1;
residual(3)= lhs-rhs;
lhs =exp(y(10));
rhs =1;
residual(4)= lhs-rhs;
lhs =exp(y(16));
rhs =T61;
residual(5)= lhs-rhs;
lhs =exp(y(17));
rhs =1-params(6)+exp(y(17))*params(6)*params(1)*exp(y(10))*exp(y(19));
residual(6)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(17))/(params(29)-exp(y(16)));
residual(7)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(12))+(exp(y(11))-exp(y(12)))*exp(y(18));
residual(8)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(19));
residual(9)= lhs-rhs;
lhs =exp(y(8))*exp(y(3));
rhs =exp(y(18))*exp(y(13));
residual(10)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(14))+exp(y(15));
residual(11)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(13))*params(6)*exp(y(19))*exp((-x(4)));
residual(12)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(3))*exp(y(8))*params(28)*exp(y(31));
residual(13)= lhs-rhs;
lhs =exp(y(11));
rhs =T119/exp(y(8));
residual(14)= lhs-rhs;
lhs =exp(y(2));
rhs =T129*exp(y(4))^(1-params(7));
residual(15)= lhs-rhs;
lhs =exp(y(8));
rhs =1;
residual(16)= lhs-rhs;
lhs =exp(y(28));
rhs =params(32)+params(31)/(1+params(5))*exp(y(22))^(1+params(5));
residual(17)= lhs-rhs;
lhs =T149;
rhs =T153;
residual(18)= lhs-rhs;
lhs =y(29);
rhs =exp(y(5))-exp(y(3))*exp(y(31))*exp(y(28));
residual(19)= lhs-rhs;
lhs =exp(y(3));
rhs =y(29)+exp(y(3))*exp(y(31));
residual(20)= lhs-rhs;
lhs =exp(y(7));
rhs =params(33)*exp(y(32));
residual(21)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(5))+exp(y(6))+exp(y(7));
residual(22)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(1))*exp(y(24));
residual(23)= lhs-rhs;
lhs =exp(y(24));
rhs =T189+(1-params(11))*T198^((-params(10))/(1-params(11)));
residual(24)= lhs-rhs;
lhs =exp(y(23));
rhs =1/exp(y(21));
residual(25)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(21))*exp(y(1))+T218;
residual(26)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(1))+T230;
residual(27)= lhs-rhs;
lhs =exp(y(34));
rhs =T238;
residual(28)= lhs-rhs;
lhs =T181^(1-params(10));
rhs =params(11)*T228+(1-params(11))*exp(y(34))^(1-params(10));
residual(29)= lhs-rhs;
lhs =exp(y(27));
rhs =exp(y(12))*T181;
residual(30)= lhs-rhs;
lhs =exp(y(27));
rhs =T251*T261*exp(x(5));
residual(31)= lhs-rhs;
lhs =y(30);
rhs =y(30)*params(18)-x(1);
residual(32)= lhs-rhs;
lhs =y(31);
rhs =y(31)*params(16)-x(2);
residual(33)= lhs-rhs;
lhs =y(32);
rhs =y(32)*params(20)-x(3);
residual(34)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(34, 34);

  %
  % Jacobian matrix
  %

  g1(1,1)=T18;
  g1(1,4)=(-(exp(y(21))*(1-params(7))*exp(y(1))*exp(y(4))))/(exp(y(4))*exp(y(4)))-(1-params(35))*(-(exp(y(6))*(-exp(y(4)))))/((1-exp(y(4)))*(1-exp(y(4))))/params(35);
  g1(1,6)=(-T26);
  g1(1,21)=T18;
  g1(2,4)=(-(params(35)*exp(y(6))^(params(35)-1)*(-exp(y(4)))*getPowerDeriv(1-exp(y(4)),1-params(35),1)));
  g1(2,6)=(-((1-exp(y(4)))^(1-params(35))*params(35)*exp(y(6))*getPowerDeriv(exp(y(6)),params(35)-1,1)));
  g1(2,9)=exp(y(9));
  g1(3,10)=params(1)*exp(y(12))*exp(y(10));
  g1(3,12)=params(1)*exp(y(12))*exp(y(10));
  g1(4,10)=exp(y(10));
  g1(5,10)=(-T61);
  g1(5,11)=(-(exp(y(10))*params(1)*(1-params(6))*exp(y(11))));
  g1(5,12)=(-(exp(y(10))*params(1)*(1-params(6))*(-exp(y(12)))));
  g1(5,16)=exp(y(16))-exp(y(16))*params(6)*params(1)*exp(y(10))*exp(y(20));
  g1(5,20)=(-(exp(y(16))*params(6)*params(1)*exp(y(10))*exp(y(20))));
  g1(6,10)=(-(exp(y(17))*params(6)*params(1)*exp(y(10))*exp(y(19))));
  g1(6,17)=exp(y(17))-exp(y(17))*params(6)*params(1)*exp(y(10))*exp(y(19));
  g1(6,19)=(-(exp(y(17))*params(6)*params(1)*exp(y(10))*exp(y(19))));
  g1(7,16)=(-((-(exp(y(17))*(-exp(y(16)))))/((params(29)-exp(y(16)))*(params(29)-exp(y(16))))));
  g1(7,17)=(-(exp(y(17))/(params(29)-exp(y(16)))));
  g1(7,18)=exp(y(18));
  g1(8,11)=(-(exp(y(11))*exp(y(18))));
  g1(8,12)=(-(exp(y(12))+exp(y(18))*(-exp(y(12)))));
  g1(8,18)=(-((exp(y(11))-exp(y(12)))*exp(y(18))));
  g1(8,19)=exp(y(19));
  g1(9,19)=(-exp(y(19)));
  g1(9,20)=exp(y(20));
  g1(10,3)=exp(y(8))*exp(y(3));
  g1(10,8)=exp(y(8))*exp(y(3));
  g1(10,13)=(-(exp(y(18))*exp(y(13))));
  g1(10,18)=(-(exp(y(18))*exp(y(13))));
  g1(11,13)=exp(y(13));
  g1(11,14)=(-exp(y(14)));
  g1(11,15)=(-exp(y(15)));
  g1(12,13)=(-(exp(y(13))*params(6)*exp(y(19))*exp((-x(4)))));
  g1(12,14)=exp(y(14));
  g1(12,19)=(-(exp(y(13))*params(6)*exp(y(19))*exp((-x(4)))));
  g1(13,3)=(-(exp(y(3))*exp(y(8))*params(28)*exp(y(31))));
  g1(13,8)=(-(exp(y(3))*exp(y(8))*params(28)*exp(y(31))));
  g1(13,15)=exp(y(15));
  g1(13,31)=(-(exp(y(3))*exp(y(8))*params(28)*exp(y(31))));
  g1(14,2)=T286;
  g1(14,3)=(-((-(exp(y(3))*exp(y(21))*params(7)*exp(y(2))))/(exp(y(3))*exp(y(3)))/exp(y(8))));
  g1(14,8)=(-((exp(y(8))*exp(y(8))*exp(y(31))-exp(y(8))*T119)/(exp(y(8))*exp(y(8)))));
  g1(14,11)=exp(y(11));
  g1(14,21)=T286;
  g1(14,28)=(-(exp(y(31))*(-exp(y(28)))/exp(y(8))));
  g1(14,31)=(-(exp(y(31))*(exp(y(8))-exp(y(28)))/exp(y(8))));
  g1(15,2)=exp(y(2));
  g1(15,3)=T298;
  g1(15,4)=(-(T129*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(7),1)));
  g1(15,22)=T298;
  g1(15,30)=(-(T129*exp(y(4))^(1-params(7))));
  g1(15,31)=T298;
  g1(16,8)=exp(y(8));
  g1(17,22)=(-(params(31)/(1+params(5))*exp(y(22))*getPowerDeriv(exp(y(22)),1+params(5),1)));
  g1(17,28)=exp(y(28));
  g1(18,2)=T149;
  g1(18,3)=(-T153);
  g1(18,21)=T149;
  g1(18,22)=(-(exp(y(21))*params(7)*exp(y(2))*exp(y(22))))/(exp(y(22))*exp(y(22)))-exp(y(3))*exp(y(31))*params(31)*exp(y(22))*getPowerDeriv(exp(y(22)),params(5),1);
  g1(18,31)=(-T153);
  g1(19,3)=exp(y(3))*exp(y(31))*exp(y(28));
  g1(19,5)=(-exp(y(5)));
  g1(19,28)=exp(y(3))*exp(y(31))*exp(y(28));
  g1(19,29)=1;
  g1(19,31)=exp(y(3))*exp(y(31))*exp(y(28));
  g1(20,3)=exp(y(3))-exp(y(3))*exp(y(31));
  g1(20,29)=(-1);
  g1(20,31)=(-(exp(y(3))*exp(y(31))));
  g1(21,7)=exp(y(7));
  g1(21,32)=(-(params(33)*exp(y(32))));
  g1(22,1)=exp(y(1));
  g1(22,5)=(-exp(y(5)));
  g1(22,6)=(-exp(y(6)));
  g1(22,7)=(-exp(y(7)));
  g1(23,1)=(-(exp(y(1))*exp(y(24))));
  g1(23,2)=exp(y(2));
  g1(23,24)=(-(exp(y(1))*exp(y(24))));
  g1(24,24)=exp(y(24))-T189;
  g1(24,33)=(-(T188*exp(y(24))*params(11)*T181*getPowerDeriv(T181,(-params(12))*params(10),1)+T187*T181*getPowerDeriv(T181,params(10),1)+(1-params(11))*(-(T181^(params(11)-1)*params(11)*T181*getPowerDeriv(T181,params(12)*(1-params(11)),1)+params(11)*T181^(params(12)*(1-params(11)))*T181*getPowerDeriv(T181,params(11)-1,1)))/(1-params(11))*getPowerDeriv(T198,(-params(10))/(1-params(11)),1)));
  g1(25,21)=(-((-exp(y(21)))/(exp(y(21))*exp(y(21)))));
  g1(25,23)=exp(y(23));
  g1(26,1)=(-(exp(y(21))*exp(y(1))));
  g1(26,10)=(-T218);
  g1(26,21)=(-(exp(y(21))*exp(y(1))));
  g1(26,25)=exp(y(25))-T218;
  g1(26,33)=(-(exp(y(25))*(T181^(params(12)*(-params(10)))*exp(y(10))*params(1)*params(11)*T181*getPowerDeriv(T181,params(10),1)+T188*exp(y(10))*params(1)*params(11)*T181*getPowerDeriv(T181,params(12)*(-params(10)),1))));
  g1(27,1)=(-exp(y(1)));
  g1(27,10)=(-T230);
  g1(27,26)=exp(y(26))-T230;
  g1(27,33)=(-(exp(y(26))*(T228*exp(y(10))*params(1)*params(11)*T181*getPowerDeriv(T181,params(10)-1,1)+exp(y(10))*params(1)*params(11)*T181^(params(10)-1)*T181*getPowerDeriv(T181,params(12)*(1-params(10)),1))));
  g1(28,25)=(-T238);
  g1(28,26)=(-(T181*(-(exp(y(26))*exp(y(25))*params(10)/(params(10)-1)))/(exp(y(26))*exp(y(26)))));
  g1(28,33)=(-T238);
  g1(28,34)=exp(y(34));
  g1(29,33)=T181*getPowerDeriv(T181,1-params(10),1)-params(11)*T181*getPowerDeriv(T181,params(12)*(1-params(10)),1);
  g1(29,34)=(-((1-params(11))*exp(y(34))*getPowerDeriv(exp(y(34)),1-params(10),1)));
  g1(30,12)=(-(exp(y(12))*T181));
  g1(30,27)=exp(y(27));
  g1(30,33)=(-(exp(y(12))*T181));
  g1(31,23)=(-(exp(x(5))*T251*T255*T256*getPowerDeriv(T256,params(14),1)*T390));
  g1(31,27)=exp(y(27))-exp(x(5))*T261*exp(y(27))*getPowerDeriv(exp(y(27)),params(15),1);
  g1(31,33)=(-(exp(x(5))*T251*T390*T256^params(14)*1/params(1)*T181*getPowerDeriv(T181,params(13),1)));
  g1(32,30)=1-params(18);
  g1(33,31)=1-params(16);
  g1(34,32)=1-params(20);
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
