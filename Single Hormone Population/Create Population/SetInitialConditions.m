function xm=SetInitialConditions(glucoseSetPoint,Num_States,modelParams)
% This function runs the plant into a glucose set-point (e.g. 160 mg/dL)
% Inputs
% Set_Point_Init:  Glucose value that states should run in to
% Num_States: Number of states in the model
% Outputs:
% xm: vector of states
Set_Point_init = glucoseSetPoint/18*0.16; 
XX0 = 0*ones(Num_States,1); 

%Q1 = Set_Point_init;

fun = @(XX) SolveSS(XX, Set_Point_init, modelParams);
xm = fsolve(fun,XX0); 

end