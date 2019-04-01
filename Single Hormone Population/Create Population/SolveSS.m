function F = SolveSS(X, Q1, ModPar)
% This function implements the differential equations for the 
% glucoregulatory model.
% Copyright 2019
% Free to use for academic research and teaching applications
% For commercial use, please contact the authors for licensing.
%
% Authors: Navid Resalat, Peter Jacobs
% Last updated: 1-30-2019
%
% Please reference using the following: 
% A statistical virtual patient population for the
% glucoregulatory system in type 1 diabetes with integrated exercise model
% Resalat N, El Youusef J, Tyler N, Castle J, Jacobs PG
% 2019
% Ui = X(1); Q2 = X(2); S1 = X(3); S2 = X(4); I = X(5); x1 = X(6); x2 = X(7); x3 = X(8);

F(1) = -X(6)*Q1 - ModPar(1) + ModPar(3)*X(2) + ModPar(6)*(1-X(8));
F(2) = X(6)*Q1 - ModPar(3)*X(2) - X(7)*X(2);

F(3) = -ModPar(10)*X(6) + (ModPar(13)*ModPar(10))*X(5);
F(4) = -ModPar(11)*X(7) + (ModPar(14)*ModPar(11))*X(5);
F(5) = -ModPar(12)*X(8) + (ModPar(15)*ModPar(12))*X(5);

F(6) = X(1) - X(3)/ModPar(7);
F(7) = X(3)/ModPar(7) - X(4)/ModPar(7);
F(8) = X(4)/(ModPar(7)*ModPar(9)) - ModPar(8) * X(5);

end