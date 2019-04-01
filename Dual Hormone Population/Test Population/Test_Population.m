% This is a test file showing the glucose data for a sample meal sceanrio and
% a sample virtual patient using the equations reporterd in Resalat et al.
% 2019, "A statistical virtual patient population for the glucoregulatory system
% in type 1 diabetes with integrated exercise model" paper, for the "Dual-
% HORMONE" control analysis. The glucose data shown here is computed with
% open-loop basal insulin delivery. you can replace the insuin input with
% any controllers' (MPC vs PID vs Fuzzy logic) output for AP in-silico analyses.
% Note: the states of the single hormone system are  [Q1;Q2;S1;S2;I;x1;x2;x3;x1g;x2g;x3g;Y;Z]
clear, clc, close all
%% loading the VPP
load('DualHormone_VPP.mat');
Ts = 5;                          % Ts: sampling interval (min);
Days_Sim = 1;                    % Number of simulation Days
Sim_time = Days_Sim*1440/Ts;     % simulation time (sample)
nn = 2;                          % ID of the virtual subject

%% loading the model parameters of the selected virtual patient

ModPar(1) = YmSbjt(1,nn);     % Fc01: Non-insulin mediated glucose uptake above 4.5 mmol/L ([mmol/kg]/min);
ModPar(2) = YmSbjt(2,nn);     % Vdg: Volume of distribution of glucose (L/kg)
ModPar(3) = YmSbjt(3,nn);     % k12: Rate constant for glucose transfer from Q2 to Q1 (min^-1)
ModPar(4) = YmSbjt(4,nn);     % Ag: Carb bioavailability (unitless);
ModPar(5) = YmSbjt(5,nn);     % tmaxG: time-to-maximum of carb absorption (min);
ModPar(6) = YmSbjt(6,nn);     % EGP0: Endogenous glucose production maximum ([mmol/kg]/min)
ModPar(7) = YmSbjt(7,nn);     % tmaxI: time-to-maximum of rapid-acting insulin absorption
ModPar(8) = YmSbjt(8,nn);     % Ke: Elimination rate of insulin (min^-1)
ModPar(9) = YmSbjt(9,nn);     % VdI: Volume of distribution of insulin (L/kg)
ModPar(10) = YmSbjt(10,nn);   % ka1: Rate constant for elimination of insulin effect from x1 (min^-1)
ModPar(11) = YmSbjt(11,nn);   % ka2: Rate constant for elimination of insulin effect from x2 (min^-1)
ModPar(12) = YmSbjt(12,nn);   % ka3: Rate constant for elimination of insulin effect from x3 (min^-1)
ModPar(13) = YmSbjt(13,nn);   % Sf1: Sensitivity factor for glucose distribution (x1) ([mU.L.min]^-2)
ModPar(14) = YmSbjt(14,nn);   % Sf2: Sensitivity factor for insulin mediated glucose utilization (x2) ([mU.L.min]^-2)
ModPar(15) = YmSbjt(15,nn);   % Sf3: Sensitivity factor for suppression of endogenous glucose production (x3) ([mU.L.min]^-1)
ModPar(16) = YmSbjt(16,nn);   % kg1: constant transfer rates [min-1]
ModPar(17) = YmSbjt(17,nn);   % kge1: elimination rates of glucagon from inaccessible [min-1]
ModPar(18) = YmSbjt(18,nn);   % kge2: elimination rates of glucagon from accessible [min-1]
ModPar(19) = YmSbjt(19,nn);   % VdGG: glucagon volume of distribution [L/kg]
ModPar(20) = YmSbjt(20,nn);   % SfGG: glucagon sensitivity factor [(ng/L)-1.min-1]
ModPar(21) = YmSbjt(21,nn);   % kc:  clearance rate of glucagon from the remote compartment [min-1]
ModPar(22) = YmSbjt(22,nn);   % kg3: constant rates [min]
ModPar(23) = YmSbjt(23,nn);   % kg2: constant transfer rates [min-1]
ModPar(24) = Weights(1,nn);   % Weight of the virtual subject (kg);
Weight = ModPar(24);
%% Select a sample meal scenario

Scenario(:,1) = [ 8; 50; 121; 303; 397; 404; 433; 566; 645; 703; 871; 914; 985];     % Time of meal event (sample)
Scenario(:,2) = [35; 79; 117; 40;  15;  100; 30;  100; 100; 100; 35;  79;  117];     % Amount of meal event (gram)
Time_Meal = round(Scenario(:,1)); Amt_Meal = Scenario(:,2);
Meal_Vector = zeros(1,Sim_time); Meal_Vector(Time_Meal) = Amt_Meal;
meal_time= []; meal_Amount = [];

%% Bolus calculations
ICR = (1700/TDIRlist(1,nn)/3);                    % ICR: Insulin to Carb Ratio
Ip = 1;                                         % percentage of pre-meal bolus  [unitless: 0-1]
Bolus = Ip*(Meal_Vector/ICR)*1000/(Weight*Ts);    % pre-meal bolus insulin        [mU/kg/min]
%% Initial Condition; the steady state run for the initial conditions
CGM_Start = 160;         % starting glucose (mg/dl)
Vdg = ModPar(2); Q1 = CGM_Start/18*Vdg;  % conversion to (mmol/kg); (mg/dl) divided by 18 (mmol/L) times Vdg (L/kg) = mmol/kg
fun = @(XX) SolveSS(XX, Q1, ModPar);   % Solve the steady states
Num_States_Plant = 13;
XX0 = zeros(8,1); XX = fsolve(fun,XX0);
xm_Plnt = [Q1; XX(2:end); zeros(5,1)];

%% Exercise Information
PVO2max = 60;  % PVO2max: Percentage of maximum oxygen consumption [0-100]
PAMM = 0.5;    % PAMM: Percentage of active muscular mass          [0-1]
Ex_Onset_time = 1*60 / Ts;    % n = 1 hours after the start of the simulation  
Ex_duration = 45 / Ts;        % Duration of the exericse bout (e.g. 45 minutes)
PGUA_1_Act = 0;               % Initial Periphery Glucose Uptake at active tissues 


%% Params of Rescue carb Algorithm
Thr_Resc = 70; Carbs_resc = 20; Win_Resc = 20/Ts;
IIR_Red_Resc = 0.25; Timer_Resc = 40/Ts;
tmax_resc = 20;
% Once Glucose < Thr_Resc [mg/dl], Carbs_resc is given [gram] with a Delay_Resc [min], ...
% glucose is remeasured after Win_Resc [min]; if Glucose is still < Thr_Resc
% another Carbs_resc will be given and glucose is remeasured after Win_Resc until Glucose > Thr_Resc;
% meanwhile Insulin insulin rates are turned down to IIR_Red_Resc percent for Timer_Resc minutes
% Note: time-to-maximum rescure carb absorption is tmax_resc [min]

Cntr_Resc = Win_Resc;
% Cntr_Resc [min]: a control flag to identify the onset of the next glucose
% measurement after a rescue carb is given; should be set to Win_Resc
resc_trig_Cntr = 0; time_resc = []; Ur_Plnt = 0; Ur_Mdl = 0;
Ins_Adj_Resc = ones(1,Sim_time+1);
%%
Ins_Basal = (TDIRlist(1,nn)/TDIR_Basal_Rate/24)*1000/Weight/60; % basal insulin (mU/kg/min)
Glu_Basal = 0;                                    % Glucagon rate [mg/kg/min]
for kk = 0:Sim_time
    
    %% Delivery of glucagon when Glucose < 90 (a simple example of glucagon delivery)
    if xm_Plnt(1,1) < (90/18*0.16)
        Glu_Basal = 0.02/Weight/60;
    else
        Glu_Basal = 0;
    end
    %%
    if kk == 0
        Ins_Total = Ins_Basal;
        Glu_Total = Glu_Basal;
    else
        Ins_Total = Ins_Adj_Resc(1,kk)*Ins_Basal + Bolus(1,kk);
        Glu_Total = Glu_Basal;
    end
    %% Meal Response
    MealResponse
    %% Exercise Response
    [M_E_PIU, M_E_PGU, M_E_HGP, PGUA_1_Act] = GetExerciseResponse(kk, PVO2max, PAMM, Ex_Onset_time, Ex_duration, PGUA_1_Act);
    
    
    %% Glucoregulatory Model
    [Ap,Bp,Cp,Dp] = DH_Glucoregulatory(xm_Plnt, ModPar, Ts, M_E_PIU, M_E_PGU, M_E_HGP);
    xm_Plnt = Ap*xm_Plnt + Bp*[Ins_Total; Glu_Total] + Dp + Ml_Vec_Plnt*(Ug_Plnt + Ur_Plnt)*Ts;    % y: mmol/kg
    Y_Plant = ((Cp*xm_Plnt)/Vdg)*18;
    BG_Output(kk+1) = Y_Plant;
    Ins_input(kk+1) = Ins_Total * Weight * 60 / 1000;  % convert from mu/kg/min to u/hr
    Glu_input(kk+1) = Glu_Total * Weight * 60 * 1000;  % conversion to mcg/hr
    Rescue_Carb_Alg;
end
Font_size = 14;
Xaxis_time = 0:(length(BG_Output)-1); Xaxis_time = Xaxis_time * Ts / 60 /24;
figure;
subplot(311); plot(Xaxis_time, BG_Output);
xhandle = xlabel('time [Day]'); yhandle = ylabel('Glucose [mg/dl]');
set(xhandle,'Fontsize',Font_size) ; set(xhandle,'Fontname','Timesnewroman');
set(yhandle,'Fontsize',Font_size) ; set(yhandle,'Fontname','Timesnewroman');
axis([0 Xaxis_time(1,end) 50 max(BG_Output)+50])
box('off')
subplot(312); plot(Xaxis_time, Ins_input);
xhandle = xlabel('time [Day]'); yhandle = ylabel('Insulin [u/hr]');
set(xhandle,'Fontsize',Font_size) ; set(xhandle,'Fontname','Timesnewroman');
set(yhandle,'Fontsize',Font_size) ; set(yhandle,'Fontname','Timesnewroman');

subplot(313); plot(Xaxis_time, Glu_input);
xhandle = xlabel('time [Day]'); yhandle = ylabel('Glucagon [mcg/hr]');
set(xhandle,'Fontsize',Font_size) ; set(xhandle,'Fontname','Timesnewroman');
set(yhandle,'Fontsize',Font_size) ; set(yhandle,'Fontname','Timesnewroman');
