%%
clear;
clc
ss.Y =log(0.84878550);
ss.Ym =log(0.84878550);
ss.K = log(5.66157077);
ss.Keff=log(5.66157077);               
ss.L =log(0.33333333);              
ss.I =log(0.14153927);             
ss.C =log(0.53748913);               
ss.G =log(0.16975710);               
ss.Q =log(1.00000000);           
ss.Lambda =log(1.00000000);            
ss.Rk =log(1.01260101);            
ss.R =log(1.01010101);               
ss.N =log(1.41539269);             
ss.Ne =log(1.40277996);             
ss.Nn =log(0.01261274);           
ss.nu =log(0.00373978);            
ss.eta =log(1.51102084);           
ss.phi =log(4.00000000);           
ss.z=log(1.02010101);              
ss.x =log(1.02010101);               
ss.Pm =log(0.76001920);                        
ss.w =log(1.29663748);            
ss.VMPK =log(0.03760101);        
ss.U =log(1.00000000);             
ss.X =log(1.31575624);            
ss.D =log(1.00000000);             
ss.F =log(2.81958684);               
ss.Z =log(3.70988896);              
ss.i =log(1.01010101);               
ss.prem =log(1.00247500);            
ss.delta =log(0.02500000);            
ss.In =0.00000000;                       
ss.Welf =-296.62068821;               
ss.a =0.00000000;               
ss.ksi =0.00000000;               
ss.g =0.00000000;               
ss.infl =0.00000000;               
ss.inflstar =0.00000000;
%ss.varrho = log(1.94246544);

variables = fieldnames(ss); 

values = 2;
for SName = variables'
    values(end+1) = ss.(SName{1});
end
values = values(2:end)'; 

variablesSS = setdiff(variables ,{'VMPK', 'Keff', 'w', 'Welf', 'prem', 'Lambda', 'R', 'X', 'a', 'ksi', 'g', 'z'});
steady_paper= 2;
for SName = variablesSS'
   if (strmatch(SName, 'In', 'exact'))
        steady_paper(end+1) = ss.(SName{1});
   else
    steady_paper(end+1) = exp(ss.(SName{1}));
   end
end
steady_paper = steady_paper(2:end)'; 
initial = steady_paper;
tab = table(variablesSS, initial);
openvar('tab');
%%
[SS_csolve, rc]=csolve(@GK_ML_stst_csolve,initial,[],1e-12,1000);
tab = table(variablesSS,SS_csolve, valueCSolve(SS_csolve));
openvar('tab');
% [SS_csolve2, rc]=csolve(@GK_ML_stst_csolve,SS_csolve,[],1e-12,10000);
% valueCSolve(SS_csolve2)
%%
tab = table(variablesSS,SS_csolve);
openvar('tab');
%%
[SS_fsolve,fval,exitflag]=fsolve(@FA_EZ_stst_fsolve,SS_csolve);

%%
diff = SS_csolve - SS_fsolve;
outTable = table(variablesSS,  SS_fsolve, SS_csolve, diff);
openvar('outTable');