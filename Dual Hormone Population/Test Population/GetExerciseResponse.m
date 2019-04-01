function [M_E_PIU, M_E_PGU, M_E_HGP, PGUA_1_Act] = GetExerciseResponse(kk, PVO2max, PAMM, Ex_Onset_time, Ex_duration, PGUA_1_Act)
strt_Ex_1 = Ex_Onset_time; End_Ex_1 = strt_Ex_1 + Ex_duration; 
Mean_PGUA_1 = 0.006*(PVO2max)^2 + 1.2264*(PVO2max) - 10.1958;

t = kk; Ts = 5;
if t >= strt_Ex_1 && t <= End_Ex_1 %&& Cntr_Ex == 1
    syms yy_1(kk); YY_1 = dsolve(diff(yy_1) == -(1/30)*yy_1 + (1/30)*Mean_PGUA_1, yy_1(0) == 0.1);
    kk = Ts*(t - strt_Ex_1);
    PGUA_1_Act = subs(YY_1); PGUA_1_Act = double(PGUA_1_Act);
    M_E_PGU = 1 + PGUA_1_Act*PAMM/35;
    M_E_PIU = 1 + 2.4*PAMM;
    M_E_HGP = 1 + PGUA_1_Act*PAMM/155;
elseif t > End_Ex_1 && t < End_Ex_1 + 5*60/Ts %&& Cntr_Ex == 1
    syms zz_1(kk); ZZ_1 = dsolve(diff(zz_1) == -(1/30)*zz_1 , zz_1(0) == PGUA_1_Act);
    kk = Ts*(t - End_Ex_1);
    PGUA_1_NoAct = subs(ZZ_1); PGUA_1_NoAct = double(PGUA_1_NoAct);
    M_E_PGU = 1 + PGUA_1_NoAct*PAMM/35;
    M_E_PIU = 1;
    M_E_HGP = 1 + PGUA_1_NoAct*PAMM/155;
else
%     Cntr_Ex = 0;
    M_E_PGU = 1;
    M_E_PIU = 1;
    M_E_HGP = 1;
end
