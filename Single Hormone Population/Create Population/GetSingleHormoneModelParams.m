function [ModPar] = GetSingleHormoneModelParams(weight)
% This function defines the model parameters for a single-hormone 
% glucoregulatory model.  The glucose kinetics and insulin kinetics are
% taken from Hovorka et al. 2004
% The insulin dynamics are taken from Wilinska et al. 2005
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
%
% For most recent code updates, please visit the AIMS lab web site
% www.ohsu.edu/jacobs
%
%% Glucose Kinetics
Fc01 = 0.0097; VdG = 0.16; k12 = 0.066;
Ag = 0.8; tmaxG = 40; EGP0 = 0.0161; 
SHparams(1,1) = Fc01; SHparams(2,1) = VdG; SHparams(3,1) = k12;
SHparams(4,1) = Ag; SHparams(5,1) = tmaxG; SHparams(6,1) = EGP0;
%% Insulin Kinetic from Hovorka et al. 2004
tmaxI = 55; Ke = 0.138; VdI = 0.12;
SHparams(7,1) = tmaxI; SHparams(8,1) = Ke; SHparams(9,1) = VdI;
%% Insulin Dynamic
ka1 = 0.006; ka2 = 0.06; ka3 = 0.03;
Sf1 = 51.2e-4; Sf2 = 8.2e-4; Sf3 = 520e-4;
SHparams(10,1) = ka1; SHparams(11,1) = ka2; SHparams(12,1) = ka3; 
SHparams(13,1) = Sf1; SHparams(14,1) = Sf2; SHparams(15,1) = Sf3;
SHparams(16,1) = weight;

ModPar = SHparams;

end