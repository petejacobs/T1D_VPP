function [Ur_Plnt, Ur_Mdl, Dr_Plnt,Dr_Mdl, Cntr_Resc,time_resc,resc_trig_Cntr] = GetRescueCarb(Y_Plant,Ur_Plnt,Ur_Mdl,Dr_Plnt,Dr_Mdl,indx,Ts,ModPar,Cntr_Resc,time_resc,resc_trig_Cntr,Thr_Resc,Carbs_resc,Win_Resc,IIR_Red_Resc,Timer_Resc,Weight,tmax_resc,delay_rescue_val)
% This function implements an adjustment to the states following a
% hypoglycemic event whereby rescue carbohydrates are given to the patient.
% This function returns the amount of updated carbohydrates taken by the
% patient as well as handles the model for absorption of these
% carbohydrates in the gut into the blood.
% Inputs:
%
% Outputs:
Ag = ModPar(4);  
if Y_Plant < (Thr_Resc)  && Cntr_Resc == Win_Resc
    
    if ~isempty(time_resc) && indx - time_resc(1,end) > 60/Ts
        Delay_Resc = delay_rescue_val/Ts;
        Ins_Adj_Resc(1,indx + Delay_Resc:indx + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc indx  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        %Resc_Carbs_subs(nn, indx + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    elseif ~isempty(time_resc) && indx - time_resc(1,end) <= Win_Resc
        Delay_Resc = 0/Ts;
        Ins_Adj_Resc(1,indx + Delay_Resc:indx + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc indx  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        %Resc_Carbs_subs(nn, indx + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    end
    
    if isempty(time_resc)
        Delay_Resc = delay_rescue_val/Ts;
        Ins_Adj_Resc(1,indx + Delay_Resc:indx + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc indx  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        %Resc_Carbs_subs(nn, indx + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    end
    
end

if resc_trig_Cntr == 1
    Ur_Plnt = (Dr_Plnt*Ag*Ts*(indx - time_resc + 1).*exp(-Ts*(indx - time_resc + 1)/tmax_resc))/(tmax_resc^2);
    Ur_Mdl = (Dr_Mdl*Ag*Ts*(indx - time_resc + 1).*exp(-Ts*(indx - time_resc + 1)/tmax_resc))/(tmax_resc^2);
    Ur_Plnt(Ur_Plnt < 0 ) = 0; Ur_Mdl(Ur_Mdl < 0 ) = 0;
end

Ur_Plnt = sum(Ur_Plnt); Ur_Mdl = sum(Ur_Mdl);

if indx > 0
    if isvector(time_resc) && (indx - time_resc(end)) >= Win_Resc - 1
        Cntr_Resc = Win_Resc;
    end
end

