% This is a test file showing the glucose data for a sample meal sceanrio and
% a sample virtual patient using the equations reporterd in Resalat et al.
% 2019, "A statistical virtual patient population for the glucoregulatory system 
% in type 1 diabetes with integrated exercise model" paper, for the "SINGLE-
% HORMONE" control analysis. The glucose data shown here is computed with 
% open-loop basal insulin delivery. you can replace the insuin input 
% with any controllers' (MPC vs PID vs Fuzzy logic) output for AP in-silico analyses.
% Note: the states of the single hormone system are  [Q1; Q2; S1; S2; I; x1; x2; x3]

clear, clc, %close all
%% loading the VPP
load('SingleHormone_VPP.mat'); 
addpath('..\Create Population');

Ts = 5;                          % Ts: sampling interval every 5 minutes (min);
Days_Sim = 4;                    % Number of simulation Days
Sim_time = Days_Sim*1440/Ts;     % simulation time (sample) 
nn = 2;                          % ID of the virtual subject

%% loading the model parameters of the selected virtual patient 

ModPar(1) = YmSbjt(1,nn);        % Fc01: Non-insulin mediated glucose uptake above 4.5 mmol/L ([mmol/kg]/min);
ModPar(2) = YmSbjt(2,nn);        % Vdg: Volume of distribution of glucose (L/kg)
ModPar(3) = YmSbjt(3,nn);        % k12: Rate constant for glucose transfer from Q2 to Q1 (min^-1)
ModPar(4) = YmSbjt(4,nn);        % Ag: Carb bioavailability (unitless);
ModPar(5) = YmSbjt(5,nn);        % tmaxG: time-to-maximum of carb absorption (min);
ModPar(6) = YmSbjt(6,nn);        % EGP0: Endogenous glucose production maximum ([mmol/kg]/min)
ModPar(7) = YmSbjt(7,nn);        % tmaxI: time-to-maximum of rapid-acting insulin absorption        
ModPar(8) = YmSbjt(8,nn);        % Ke: Elimination rate of insulin (min^-1)    
ModPar(9) = YmSbjt(9,nn);        % VdI: Volume of distribution of insulin (L/kg)  
ModPar(10) = YmSbjt(10,nn);      % ka1: Rate constant for elimination of insulin effect from x1 (min^-1) 
ModPar(11) = YmSbjt(11,nn);      % ka2: Rate constant for elimination of insulin effect from x2 (min^-1)
ModPar(12) = YmSbjt(12,nn);      % ka3: Rate constant for elimination of insulin effect from x3 (min^-1) 
ModPar(13) = YmSbjt(13,nn);      % Sf1: Sensitivity factor for glucose distribution (x1) ([mU.L.min]^-2)
ModPar(14) = YmSbjt(14,nn);      % Sf2: Sensitivity factor for insulin mediated glucose utilization (x2) ([mU.L.min]^-2) 
ModPar(15) = YmSbjt(15,nn);      % Sf3: Sensitivity factor for suppression of endogenous glucose production (x3) ([mU.L.min]^-1)
ModPar(16) = Weights(1,nn);      % Weight of the virtual subject (kg); 
Weight = ModPar(16);
%% Select a sample meal scenario

Scenario(:,1) = [ 8; 50; 121; 303; 397; 404; 433; 566; 645; 703; 871; 914; 985];     % Time of meal event (sample)
Scenario(:,2) = [35; 79; 117; 40;  15;  100; 30;  100; 100; 100; 35;  79;  117];     % Amount carbs in meal event (gram)
Time_Meal = round(Scenario(:,1)); Amt_Meal = Scenario(:,2);
Meal_Vector = zeros(1,Sim_time); Meal_Vector(Time_Meal) = Amt_Meal;
meal_time= []; meal_Amount = [];

%% Bolus calculations
ICR = (1700/TDIRlist(1,nn)/3);       % ICR: Insulin to Carb Ratio
Ip = 1;                                         % percentage of pre-meal bolus  [unitless: 0-1]
Bolus = Ip*(Meal_Vector/ICR)*1000/(Weight*Ts);    % pre-meal bolus insulin with correct units for model        [mU/kg/min]
%% Initial Condition; the steady state run for the initial conditions  
CGM_Start = 160;         % starting glucose (mg/dl)
Num_States_Plant=8;
xm_Plnt=SetInitialConditions(CGM_Start,Num_States_Plant,ModPar);
Q1 = CGM_Start/18*ModPar(2);  % conversion to (mmol/kg); (mg/dl) divided by 18 (mmol/L) times Vdg (L/kg) = mmol/kg
xm_Plnt = [Q1; xm_Plnt(2:end)];  % Replace Q1 compartment with the starting glucose level

%% Exercise Information
PVO2max = 60;  % PVO2max: Percentage of maximum oxygen consumption [0-100]
PAMM = 0.5;    % PAMM: Percentage of active muscular mass          [0-1]
Ex_Onset_time = 1*60 / Ts;    % n = 1 hours after the start of the simulation  
Ex_duration = 45 / Ts;        % Duration of the exericse bout (e.g. 45 minutes)
PGUA_1_Act = 0;               % Initial Periphery Glucose Uptake at active tissues 

%% Params of Rescue carb Algorithm
tmax_resc = 20; 
% Once Glucose < Thr_Resc [mg/dl], Carbs_resc is given [gram] with a Delay_Resc [min], ...
% glucose is remeasured after Win_Resc [min]; if Glucose is still < Thr_Resc
% another Carbs_resc will be given and glucose is remeasured after Win_Resc until Glucose > Thr_Resc;
% meanwhile Insulin insulin rates are turned down to IIR_Red_Resc percent for Timer_Resc minutes 
% Note: time-to-maximum rescure carb absorption is tmax_resc [min] 
Thr_Resc = 70;  % Rescue carbs given for glucose < 70 mg/dL
Carbs_resc = 20; % 20 g of carbs given when glucose <  70 mg/dL
Win_Resc = 40/Ts;  % Window for lower insulin dosed is 40 minutes after hypo
IIR_Red_Resc = 0.25; Timer_Resc = 40/Ts;  % Insulin is reduced to 25% for 40 minutes after a hypo
delay_rescue_val = 20; % Rescue carb is given 20 minutes after hypo occurs
Cntr_Resc = Win_Resc;     
% Cntr_Resc [min]: a control flag to identify the onset of the next glucose
% measurement after a rescue carb is given; should be set to Win_Resc  
resc_trig_Cntr = 0; time_resc = []; Ur_Plnt = 0; Ur_Mdl = 0;Dr_Plnt = 0; Dr_Mdl=0;
Ins_Adj_Resc = ones(1,Sim_time+1); 
%% 

% Provide patient with constant basal insulin
u_Basal = (TDIRlist(1,nn)/TDIR_Basal_Rate/24)*1000/Weight/60; % basal insulin (mU/kg/min)

for kk = 0:Sim_time
    if kk == 0
        u_Total = u_Basal;
    else
        u_Total = Ins_Adj_Resc(1,kk)*u_Basal + Bolus(1,kk);
    end
    %% Meal Response
    [Ug_Plnt,Ml_Vec_Plnt,meal_time,meal_Amount] = GetMealResponse(ModPar,Num_States_Plant,kk,Ts,meal_Amount,meal_time,Meal_Vector,Weight);
    %% Exercise Response
    [M_E_PIU, M_E_PGU, M_E_HGP, PGUA_1_Act] = GetExerciseResponse(kk, PVO2max, PAMM, Ex_Onset_time, Ex_duration, PGUA_1_Act);
    %%
    [Ap,Bp,Cp,Dp] = SH_Glucoregulatory(xm_Plnt, ModPar, Ts, M_E_PIU, M_E_PGU, M_E_HGP);
    xm_Plnt = Ap*xm_Plnt + Bp*u_Total + Dp + Ml_Vec_Plnt*(Ug_Plnt + Ur_Plnt)*Ts;    % y: mmol/kg
    Y_Plant = ((Cp*xm_Plnt)/ModPar(2))*18;
    BG_Output(kk+1) = Y_Plant;
    Ins_input(kk+1) = u_Total * Weight * 60 / 1000;  % convert from mu/kg/min to u/hr
    [Ur_Plnt, Ur_Mdl, Dr_Plnt, Dr_Mdl, Cntr_Resc,time_resc,resc_trig_Cntr] = GetRescueCarb(Y_Plant,Ur_Plnt,Ur_Mdl,Dr_Plnt,Dr_Mdl, kk,Ts,ModPar,Cntr_Resc,time_resc,resc_trig_Cntr,Thr_Resc,Carbs_resc,Win_Resc,IIR_Red_Resc,Timer_Resc,Weight,tmax_resc,delay_rescue_val);
    %Rescue_Carb_Alg;
end
Font_size = 14;
Xaxis_time = 0:(length(BG_Output)-1); Xaxis_time = Xaxis_time * Ts / 60 /24;
figure; 
subplot(211); plot(Xaxis_time, BG_Output); 
xhandle = xlabel('time [Day]'); yhandle = ylabel('Glucose [mg/dl]');
set(xhandle,'Fontsize',Font_size) ; set(xhandle,'Fontname','Timesnewroman');
set(yhandle,'Fontsize',Font_size) ; set(yhandle,'Fontname','Timesnewroman');
axis([0 Xaxis_time(1,end) 50 max(BG_Output)+50])
box('off')
subplot(212); plot(Xaxis_time, Ins_input);
xhandle = xlabel('time [Day]'); yhandle = ylabel('Insulin [u/hr]');
set(xhandle,'Fontsize',Font_size) ; set(xhandle,'Fontname','Timesnewroman');
set(yhandle,'Fontsize',Font_size) ; set(yhandle,'Fontname','Timesnewroman');
