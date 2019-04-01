% Create_DH_VPP.m
% This script generates a population of virtual patients that is 
% represented by a mathematical glucoregulatory model.  The model includes
% a representation of insulin kinetics, insulin dynamics, and carbohydrate
% kinetics.
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
%% Add the common functions folder to the directory
addpath('..\..\Common');
%% set the clinical average TDIR and weight of Type 1 diabetes
Avg_Wgt = 76.3; Ave_TDIR = 42.3; Population_No = 100; 
target_glucose = 115; % Used to set steady-state conditions
%% TDIR vs Sc - Figure 2 of the Resalat et al. 2019 paper
% This function gets the insulin kinetics, glucose kinetics and insulin 
% dynamic parameters in the model
            % DH_Params % removed this script
DHparams = GetDualHormoneModelParams(Avg_Wgt);
Weight = DHparams(end);
% Ratio between the TDIR and TDIR_basal.  
% A ratio of 2 indicates basal is 50% of TDIR (1/0.5). Ratio of 1.78 indicates 
% that basal insulin comprises 56.18% of total daily insulin (1/0.5618)
TDIR_Basal_Rate = 1.78;

% Create an array of insulin sensitivity composite (Sc) values ranging 
% from 0.05 to 2, whereby an Sc of 1.0 indicates the average insulin
% sensitivity of a person without diabetes
ModPar = DHparams;
% initialize matrices
Sc = zeros(1,40);
IIR_SS = zeros(1,40);
Init_Cond = zeros(13,40);
for Mlp = 1:40
    % set Sc in array
    Sc(:,Mlp) = Mlp*0.05; 
    
    % Using current Sc value, update the insulin sensitivity model
    % parameters (Sf1, Sf2, and Sf3)
    ModPar(13:15,:) = Sc(:,Mlp)*DHparams(13:15,1); 
    
    % Run the model, without any meals, and get the glucose to converge to a
    % target of 115 mg/dL  
    [xm_Plnt, iir_ss_temp] = GetPlantSteadyState(target_glucose,ModPar);

    % Insulin infusion rate is in mU/kg/min
    IIR_SS(:,Mlp) = iir_ss_temp(1,1)*Weight*60/1000;    % iir_ss_temp(1,1) is in mU/kg/min - conversion to U/hr
    Init_Cond(:,Mlp) = xm_Plnt;
end
TDIR = 24*TDIR_Basal_Rate*IIR_SS;   % Estimated TDIR

% Find sensitivity composite (Sc) reduction factor that is closest match
% based on TDIR determined above during steady-state
A = repmat(Ave_TDIR,[1 length(TDIR)]);
[minValue,closestIndex] = min(abs(A-TDIR));
Ins_red_fac = Sc(closestIndex); 

% Ins_red_fac shows the ratio f the TDIR of a person without diabetes to
% a person with diabetes. A ratio of 1.0 would mean that they have the 
% same insulin sensitivity.  We determined that for a TDIR of 42, the 
% Sc=0.4 yields an accurate steady-state value for total insulin consumed
% See paper Resalat et al. for more details.
%% Creating VP. based on the prior result
% Virtual_Population;  
% creates specified number of virtual patients
[HovParams, Weights] = SampleVirtualPopulationFromDistributionDH(DHparams,Ins_red_fac,Population_No,Avg_Wgt);
% Prune_Population;    
% prunes the population based on the physiological tests
[RevisedPar_4_pass,Weightsdata_4_pass]=PrunePopulationDH(HovParams,Weights);

YmSbjt = RevisedPar_4_pass; Weights = Weightsdata_4_pass;
IIR_SS = []; Init_Cond = []; Num_VP = size(YmSbjt,2);
glucoseSetPoint=115; % Set point for all virtual patients to calculate TDIR
for nn = 1:Num_VP
%     Plant_init_Cond
    xm_Plnt=SetInitialConditions(glucoseSetPoint,13,YmSbjt(:,nn));  % steady state run for the target of 115 mg/dl
    IIR_SS(:,nn) = xm_Plnt(1,1)*Weights(nn)*60/1000;    % XX(1,1) is in mU/kg/min
%     Init_Cond(:,nn) = xm_Plnt;
end
TDIRlist = 24*TDIR_Basal_Rate*IIR_SS;
figure; plot(TDIRlist,'-*'); figure; hist(TDIRlist,20)

clearvars -except TDIRlist Weights YmSbjt TDIR_Basal_Rate sub

save('DualHormone_VPP.mat')

