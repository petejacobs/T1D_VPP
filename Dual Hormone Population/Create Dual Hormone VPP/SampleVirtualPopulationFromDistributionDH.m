function [HovParams, Weights] = SampleVirtualPopulationFromDistributionDH(DHparams,Ins_red_fac,Population_No,Avg_Wgt)
% This function creates the virtual patient population by sampling the
% insulin sensitivity factor params from a distribution about the mean.
% Inputs:
% SHparams:  Model parameters for single-hormone VPP
% Ins_red_fac:  This is the amount that we are reducing the insulin
%               sensitivity relative to a person without type 1 diabetes
%               In the Resalat 2019 paper, this value was set to 0.4
% Population_No: Total number of patients to create
% Avg_Wgt:  The average weight in kg of the patient population
% Outputs:
% HovParams:  This is an array of model parameters for each of the new
%             patients created
% Weights:  This is an array of weights for each of the new patients
%           created
Parameters_No = 23; 
Sf2_Sf1_Corr = 0.75; Sf2_Sf1_STD_Ratio = 1 - Sf2_Sf1_Corr;
Sf2_Sf3_Corr = 0.25; Sf2_Sf3_STD_Ratio = 1 - Sf2_Sf3_Corr;
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')

% We describe a linear correlation between insulin sensitivity
% factors Sf1, Sf2, and Sf3 in the Resalat 2019 et al. paper
% Here we randomize the Sf2 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate according to parameters a and b.
% The insulin sensitivity reduction factor is also applied here.

a = 0.0016; b = 0.0161;
EGP0_samples = a.*randn(Population_No,1) + b;

% Randomize the Sf2 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = 0.00032*Ins_red_fac; b = 0.00082*Ins_red_fac;
rng(0,'twister'); RND_1 = randn(Population_No,1);
RND_2 = rand(Population_No,1); Sf2_samples = a.*RND_1 + b;

% Randomize the Sf1 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = 0.00131*Ins_red_fac; b = 0.00512*Ins_red_fac;
W1 = -Sf2_Sf1_STD_Ratio; W2 = Sf2_Sf1_STD_Ratio;
Sf1_samples = a.*(RND_1 + (W1 + (W2-W1).*RND_2)) + b;

% Randomize the Sf3 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = 0.0125*Ins_red_fac; b = 0.052*Ins_red_fac;
W1 = -Sf2_Sf3_STD_Ratio; W2 = Sf2_Sf3_STD_Ratio;
Sf3_samples = a.*(RND_1 + (W1 + (W2-W1).*RND_2)) + b;

row_Sf1 = find(Sf1_samples <= 0); row_Sf2 = find(Sf2_samples <= 0);
row_Sf3 = find(Sf3_samples <= 0);
row_Zeros = [row_Sf1; row_Sf2; row_Sf3]; row_Zeros = unique(row_Zeros);

Sf1_samples(row_Zeros) = []; Sf2_samples(row_Zeros) = []; Sf3_samples(row_Zeros) = [];
EGP0_samples(row_Zeros) = [];

Study_Pop = length(Sf2_samples); Params = zeros(Parameters_No,Study_Pop);

for pp = 1:Parameters_No
    Params(pp,:) =  DHparams(pp,1);
end

Params(13,:) = Sf1_samples'; Params(14,:) = Sf2_samples'; Params(15,:) = Sf3_samples';
Params(6,:) = EGP0_samples';

b = DHparams(20,1); a = b/3;
SfGG_samples = a.*randn(Population_No,1) + b;

b = DHparams(21,1); a = b/3;
kc_samples = a.*randn(Population_No,1) + b;

b = DHparams(22,1); a = b/3;
kg3_samples = a.*randn(Population_No,1) + b;

SfGG_samples(row_Zeros) = []; kc_samples(row_Zeros) = []; kg3_samples(row_Zeros) = [];

Params(20,:) = SfGG_samples'; Params(21,:) = kc_samples'; Params(22,:) = kg3_samples';
YmSbjt = Params;

a = 15.9; b = Avg_Wgt;
Wgt_samples = a.*randn(Population_No,1) + b;
Wgt_samples(row_Zeros) = [];
Weights = Wgt_samples';

col = prod(YmSbjt,1)~=0; YmSbjt = YmSbjt(:,col);
col = find(prod(Weights,1)~=0); Weights = Weights(:,col);
Weights = Weights(:,1:size(YmSbjt,2));

%% remove subs with any negative parameter
[~, col] = find(YmSbjt <= 0);
YmSbjt_temp = []; Weightsdata_temp = [];
for ll = 1:size(YmSbjt,2)
    if sum(ll == col) == 0
        YmSbjt_temp = [YmSbjt_temp YmSbjt(:,ll)];
        Weightsdata_temp = [Weightsdata_temp Weights(:,ll)];
    end
end
YmSbjt = YmSbjt_temp; Weights = Weightsdata_temp;
HovParams = YmSbjt;

end