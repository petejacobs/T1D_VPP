function [HovPars] = CreateVirtualPatientPopulation(HovParams)
Parameters_No = 15; VB = 3; Std_Ste = 1; % VB: Variation of each parameter around mean
Sf2_Sf1_Corr = 0.75; Sf2_Sf1_STD_Ratio = 1 - Sf2_Sf1_Corr;
Sf2_Sf3_Corr = 0.25; Sf2_Sf3_STD_Ratio = 1 - Sf2_Sf3_Corr;
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')

% Randomize the Sf2 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = Std_Ste*0.00032*Ins_red_fac; b = 0.00082*Ins_red_fac; 
rng(0,'twister'); RND_1 = randn(Population_No,1);
RND_2 = rand(Population_No,1);
Sf2_samples = a.*RND_1 + b; 

% Randomize the Sf1 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = Std_Ste*0.00144*Ins_red_fac; b = 0.00512*Ins_red_fac; 
W1 = -Sf2_Sf1_STD_Ratio; W2 = Sf2_Sf1_STD_Ratio;
Sf1_samples = a.*(RND_1 + (W1 + (W2-W1).*RND_2)) + b; 

% Randomize the Sf3 parameter while ensuring that Sf1, Sf2, and Sf3
% correlate
a = Std_Ste*0.01375*Ins_red_fac; b = 0.052*Ins_red_fac; 
W1 = -Sf2_Sf3_STD_Ratio; W2 = Sf2_Sf3_STD_Ratio;
Sf3_samples = a.*(RND_1 + (W1 + (W2-W1).*RND_2)) + b; 

row_Sf1 = find(Sf1_samples <= 0); row_Sf2 = find(Sf2_samples <= 0);
row_Sf3 = find(Sf3_samples <= 0);
row_Zeros = [row_Sf1; row_Sf2; row_Sf3]; row_Zeros = unique(row_Zeros);

Sf1_samples(row_Zeros) = []; Sf2_samples(row_Zeros) = []; Sf3_samples(row_Zeros) = [];
Study_Pop = length(Sf2_samples); Params = zeros(Parameters_No,Study_Pop);

for pp = 1:Parameters_No
    Params(pp,:) =  HovParams(pp,1);
end
Params(13,:) = Sf1_samples'; Params(14,:) = Sf2_samples'; Params(15,:) = Sf3_samples';
HovParams = Params;

a = 15.9; b = Avg_Wgt; Wgt_samples = a.*randn(Study_Pop,1) + b; Weights = Wgt_samples';
end
