function F = SolveSS(X, Q1, ModPar)
U_glucagon = 0;

% Ui = X(1); Q2 = X(2); 
% S1 = X(3); S2 = X(4); I = X(5);
% x1 = X(6); x2 = X(7); x3 = X(8); 
% x1g = X(9); x2g = X(10); x3g = X(11);
% Y = X(12); Z = X(13);

Fc01  = ModPar(1); VdG = ModPar(2); k12 = ModPar(3); 
Ag = ModPar(4); tmaxG =  ModPar(5); EGP0 = ModPar(6);

tmaxI = ModPar(7); Ke = ModPar(8); VdI = ModPar(9);
ka1 = ModPar(10); ka2 = ModPar(11); ka3 = ModPar(12);
Sf1 = ModPar(13); Sf2 = ModPar(14); Sf3 = ModPar(15); 

kg1 = ModPar(16); kge1 = ModPar(17); kge2 = ModPar(18);         
VdGG = ModPar(19);       
SfGG = ModPar(20); kc = ModPar(21); kg3 = ModPar(22); kg2 = ModPar(23);        
kg = ((10^6)*kc*SfGG)/VdGG;

F(1) = -X(6)*Q1 - Fc01 + k12*X(2) + EGP0*(1 - X(8));
F(2) = X(6)*Q1 - k12*X(2) - X(7)*X(2);

F(3) = X(1) - X(3)/tmaxI;
F(4) = X(3)/tmaxI - X(4)/tmaxI;
F(5) = X(4)/(tmaxI*VdI) - Ke * X(5);

F(6) = -ka1*X(6) + (Sf1*ka1)*(X(5));
F(7) = -ka2*X(7) + (Sf2*ka2)*(X(5));
F(8) = -ka3*X(8) + (Sf3*ka3)*(X(5));

% F(10) = -(kg1 + kge1)*X(10) + U_glucagon;
% F(11) = kg1*X(10) - kg2*X(11);
% F(12) = kg2*X(11) - kge2*X(12);
% 
% F(13) = kg*X(12) - kc*X(13);
% F(14) = kg*kg2*X(11) - kg*kge2*X(12) - kc*X(14); 

end