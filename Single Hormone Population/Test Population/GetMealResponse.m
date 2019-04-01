function [Ug_Plnt, Ml_Vec_Plnt,meal_time,meal_Amount] = GetMealResponse(ModPar,Num_States_Plant,indx,Ts,meal_Amount,meal_time,Meal_Vector,Weight)
% This function calculates the impact of a meal using the meal model
% described in Hovorka et al. 2004.
% Inputs
% ModPar:  Model parameters of glucoregulatory model
% indx:  Current sample index
% Ts:  sample period (minutes)
% meal_Amount:  amount of last meal consumed
% Meal_Vector:  vector of meal amounts
% Outputs
% Ug_Plnt: Amount of meal that the plant consumes
% Ml_Vec_Plnt: Vector indicating which part of state is plant
% meal_time: Time meal was consumed
% meal_amount: amoutn of last meal consumed, must be returned to calling
% function
Ag = ModPar(4); tmax_CHO = ModPar(5);
if indx > 0
    if Meal_Vector(1,indx) ~= 0
        meal_time = [meal_time indx]; meal_Amount = [meal_Amount Meal_Vector(1,indx)];
    end
end
if isempty(meal_time)
    Ug_Plnt = 0;
else
    Dg_Plnt = ((meal_Amount/Weight)/0.18);
    Ug_Plnt_Temp = (Dg_Plnt*Ag*Ts.*(indx-meal_time).*exp(-Ts*(indx-meal_time)/tmax_CHO))/(tmax_CHO^2);
    Ug_Plnt = sum(Ug_Plnt_Temp);
end
Ml_Vec_Plnt = zeros(Num_States_Plant,1); Ml_Vec_Plnt(1,1) = 1;
