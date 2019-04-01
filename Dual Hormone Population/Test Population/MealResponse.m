Ag = ModPar(4); tmax_CHO = ModPar(5);
if kk > 0
    if Meal_Vector(1,kk) ~= 0
        meal_time = [meal_time kk]; meal_Amount = [meal_Amount Meal_Vector(1,kk)];
    end
else 
    Ug_Plnt = 0;
end
if isempty(meal_time)
    Ug_Plnt = 0;
else
    Dg_Plnt = ((meal_Amount/Weight)/0.18);
    Ug_Plnt_Temp = (Dg_Plnt*Ag*Ts.*(kk-meal_time).*exp(-Ts*(kk-meal_time)/tmax_CHO))/(tmax_CHO^2);
    Ug_Plnt = sum(Ug_Plnt_Temp);
end
Ml_Vec_Plnt = zeros(Num_States_Plant,1); Ml_Vec_Plnt(1,1) = 1;
