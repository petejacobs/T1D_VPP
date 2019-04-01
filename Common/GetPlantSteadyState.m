function [plant_states_ss, IIR_SS] = GetPlantSteadyState(target_glucose,ModPar)
%
% [plant_states_ss, IIR_SS] = GetPlantSteadyState(target_glucose,ModPar)
% target_glucose: glucose in mg/dL to solve for at steady state
% ModPar: parameters of the glucoregulatory model
% plan_states_ss:  States of the model at steady-state
% IIR_SS:  insulin infusion rate in mU/kg/min
%
% This function solves for the model states and basal insulin required 
% to bring the virtual patient to a target_glucose level (e.g. 115 mg/dL)
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

VdG = ModPar(2,1);

% Convert target glucose from mg/dl to mmol/kg
Set_Point_init = target_glucose/18*VdG; 

Num_States = 13; XX0 = zeros(Num_States,1); 
fun = @(XX) SolveSS(XX, Set_Point_init, ModPar);
XX = fsolve(fun,XX0); X_Total_Plnt = XX;
xm_Plnt = [Set_Point_init; X_Total_Plnt(2:end,1)];

IIR_SS = XX;
plant_states_ss = xm_Plnt;
end