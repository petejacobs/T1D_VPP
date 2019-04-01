function [RevisedPar_2_pass, Weightsdata_2_pass, RevisedPar_2_fail,Weightsdata_2_fail]=PrunePopulation(HovParams, Weights)
% This function will prune a virtual patient population by requiring that
% each of the patients, with parameters given by HovParams passes a set of
% physiologic tests.  The tests include the following:
%
%
%
% Inputs:
% HovParams:  This lists the model parameters for each the vitual patients
% Weights: List of weight in kg of each virtual patient
% Outputs:
% RevisedPar_2_pass:  This is a list of all virtual patients who passed
% RevisedPar_2_fail:  This is a list of all virtual patients who failed
%
%% Step 1

% Initialize arrays, and indicate the time duration over which to run
% (1000 samples)
N_sim = 1000; [~, NumSub] = size(HovParams);
RevisedPar_1_pass = []; Weightsdata_1_pass = [];
RevisedPar_1_fail = []; Weightsdata_1_fail = [];
Num_fail_1 = []; Num_fail_2 = [];

% For test 1, we give the patients no insulin, and ensure that their
% glucose levels rise above 300 mg/dL
for nn = 1:NumSub
%    Sim_parameters
%    Hovorka_Parameters;
    % Select model parameters for virtual patient nn
    ModPar = HovParams(:,nn)';Weight=Weights(nn);ModPar(16)=Weight;
    
    % Initialize the states of the model
    xm = zeros(8,1);
    
    % Across N_sim time stamps, run the model presuming no insulin is given
    % to the patient
    for kk = 1:N_sim
        [Apr, Bpr, Cpr, Dpr] = Plant(xm, ModPar);
        u = 0; % this is the insulin given
        xm = Apr*xm + Bpr*u + Dpr; y = Cpr*xm; xm_sub_S1(nn,kk) = xm(1,1);
    end
    
    % We only pass this patient if their glucose is greater than 300 mg/dL
    % after giving no insulin across 1000 time steps.
    if  (xm(1,1) * 18 / ModPar(2)) >= 300  % ModPar(2) is volume of distribution (Vg)
        RevisedPar_1_pass = [RevisedPar_1_pass HovParams(:,nn)];
        Weightsdata_1_pass = [Weightsdata_1_pass Weights(1,nn)];
    else
        Num_fail_1 = [Num_fail_1 nn];
        RevisedPar_1_fail = [RevisedPar_1_fail HovParams(:,nn)];
        Weightsdata_1_fail = [Weightsdata_1_fail Weights(1,nn)];
    end
end

%% Step 2
startingGlucose=160;  % Starting glucose level to run the system in to
numOfStates = 8;  % Number of states of the glucoregulatory model
HovParams = RevisedPar_1_pass; otherdata = Weightsdata_1_pass; Weights = otherdata;
N_sim = 1000; [~, NumSub] = size(HovParams);
RevisedPar_2_pass = []; Weightsdata_2_pass = [];
RevisedPar_2_fail = []; Weightsdata_2_fail = [];

% For second test, give a large amount of insulin, and make sure glucose
% drops below 100 mg/dL
for nn = 1: NumSub
%    Sim_parameters
%    Hovorka_Parameters; 
    % Select model parameters for virtual patient nn
    ModPar = HovParams(:,nn)';Weight=Weights(nn);ModPar(16)=Weight;
    
%    Initial_Conditions
    % Initialize the states of the model to the starting glucose level
    xm = SetInitialConditions(startingGlucose,numOfStates,ModPar);
    for kk = 1:N_sim
        [Apr, Bpr, Cpr, Dpr] = Plant(xm, ModPar);
        u = 1000*(15/60/Weight);  % Dose a large amount of insulin
        xm = Apr*xm + Bpr*u + Dpr; y = Cpr*xm;
        xm(xm < 0) = 0;  xm_sub_S2(nn,kk) = xm(1,1);
    end
    
    % Make sure insulin drops below 100 mg/dL to pass
    if (xm(1,1) * 18 / ModPar(2)) <= 100  % ModPar(2) is volume of distribution (Vg)
        RevisedPar_2_pass = [RevisedPar_2_pass HovParams(:,nn)];
        Weightsdata_2_pass = [Weightsdata_2_pass otherdata(1,nn)];
    else
        Num_fail_2 = [Num_fail_2 nn];
        RevisedPar_2_fail = [RevisedPar_2_fail HovParams(:,nn)];
        Weightsdata_2_fail = [Weightsdata_2_fail otherdata(1,nn)];
    end
end
end
