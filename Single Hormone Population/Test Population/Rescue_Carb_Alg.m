Ag = ModPar(4); tmax_CHO = ModPar(5);
if Y_Plant < (Thr_Resc)  && Cntr_Resc == Win_Resc
    
    if ~isempty(time_resc) && kk - time_resc(1,end) > 60/Ts
        Delay_Resc = 10/Ts;
        Ins_Adj_Resc(1,kk + Delay_Resc:kk + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc kk  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        Resc_Carbs_subs(nn, kk + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    elseif ~isempty(time_resc) && kk - time_resc(1,end) <= Win_Resc
        Delay_Resc = 0/Ts;
        Ins_Adj_Resc(1,kk + Delay_Resc:kk + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc kk  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        Resc_Carbs_subs(nn, kk + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    end
    
    if isempty(time_resc)
        Delay_Resc = 10/Ts;
        Ins_Adj_Resc(1,kk + Delay_Resc:kk + Delay_Resc + Timer_Resc) = IIR_Red_Resc;
        resc_trig_Cntr = 1;
        time_resc = [time_resc kk  + Delay_Resc];
        Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
        Resc_Carbs_subs(nn, kk + Delay_Resc) = resc_trig_Cntr;
        Cntr_Resc = 0;
    end
    
end

if resc_trig_Cntr == 1
    Ur_Plnt = (Dr_Plnt*Ag*Ts*(kk - time_resc + 1).*exp(-Ts*(kk - time_resc + 1)/tmax_resc))/(tmax_resc^2);
    Ur_Mdl = (Dr_Mdl*Ag*Ts*(kk - time_resc + 1).*exp(-Ts*(kk - time_resc + 1)/tmax_resc))/(tmax_resc^2);
    Ur_Plnt(Ur_Plnt < 0 ) = 0; Ur_Mdl(Ur_Mdl < 0 ) = 0;
end

Ur_Plnt = sum(Ur_Plnt); Ur_Mdl = sum(Ur_Mdl);

if kk > 0
    if isvector(time_resc) && (kk - time_resc(end)) >= Win_Resc - 1
        Cntr_Resc = Win_Resc;
    end
end


% 
% if (Y_Plant < Thr_Resc) && Cntr_Resc == Win_Resc
%     
%     Ins_Adj_Resc(1,kk + 1:kk + 1 + Timer_Resc) = 0.25 * Ins_Adj_Resc(1,kk);
%     resc_trig_Cntr = 1;
%     time_resc = [time_resc kk];
%     Dr_Plnt = (Carbs_resc/Weight)/0.18; Dr_Mdl = (Carbs_resc/Weight)*1000;
%     Resc_Carbs_subs(nn, kk) = resc_trig_Cntr;
%     Cntr_Resc = 0;
%     
% end
% 
% if resc_trig_Cntr == 1;
%     Ur_Plnt = (Dr_Plnt*Ag*Ts*(kk - time_resc + 1).*exp(-Ts*(kk - time_resc + 1)/tmax_resc))/(tmax_resc^2);
%     Ur_Mdl = (Dr_Mdl*Ag*Ts*(kk - time_resc + 1).*exp(-Ts*(kk - time_resc + 1)/tmax_resc))/(tmax_resc^2);
% end
% 
% Ur_Plnt = sum(Ur_Plnt); Ur_Mdl = sum(Ur_Mdl);
% 
% if kk > 0
%     if isvector(time_resc) && (kk - time_resc(end)) >= Win_Resc
%         Cntr_Resc = Win_Resc;
%     end
% end
