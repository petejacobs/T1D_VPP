Target = 115; % (mg/dl)
Set_Point_init = Target/18*VdG; Num_States = 8; XX0 = 0*ones(Num_States,1); 
ModPar = YmSbjt(:,nn);
fun = @(XX) SolveSS(XX, Set_Point_init, ModPar);
XX = fsolve(fun,XX0); X_Total_Plnt = XX;
xm_Plnt = [Set_Point_init; X_Total_Plnt(2:end,1)];
Weight = Weights(1,nn);
