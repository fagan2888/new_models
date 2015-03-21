ss.Y =log(0.84878550);
ss.Ym =log(0.84878550);
ss.K = log(5.66157077);
ss.Keff=log(5.66157077);               
ss.L =log(0.33333333);              
ss.I =log(0.14153927);             
ss.C =log(0.53748913);               
ss.G =log(0.16975710);               
ss.Q =log(1.00000000);           
ss.varrho =log(1.94246544);         
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

variables = fieldnames(ss); 

values = 2;
for SName = variables'
    values(end+1) = ss.(SName{1});
end
values = values(2:end)'; 


%% Exponential steady as initial
initial = values;
[SS_VALUES_paper,rc]=csolve(@FA_stst_csolve_exp, initial ,[],1e-12,10000)
initial = [0.9; 0.9; values(3:end)];
[SS_VALUES_paper,rc]=csolve(@FA_stst_csolve_exp,initial,[],1e-12,10000)


%% reduced system level test steady state values
variablesSS = setdiff(variables ,{'VMPK', 'Keff', 'w', 'Welf', 'prem'});
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

[SS_opti, rc]=csolve(@FA_stst_csolve_2_trans,initial,[],1e-12,1000);
diff = SS_opti - steady_paper;
table(variablesSS, initial, SS_opti, steady_paper,diff)


%% reduced system level test steady state values
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

initial = [0.9*ones(16,1); steady_paper(17:18); 0.9* ones(3,1); steady_paper(22:23); 0.1; 1.55; 0.9; 0.9];


[SS_opti, rc]=csolve(@FA_stst_csolve_2_trans2,initial,[],1e-12,10000);
diff = SS_opti - steady_paper;
table(variablesSS, initial, SS_opti, steady_paper,diff)
max(abs(diff)) 
rc
%%
csolve(@FA_stst_csolve_2_trans2,SS_opti,[],1e-12,1);

