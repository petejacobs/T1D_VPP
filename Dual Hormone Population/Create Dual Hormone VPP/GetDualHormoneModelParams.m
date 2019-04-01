function [ModPar] = GetDualHormoneModelParams(weight)
% This function defines the model parameters for a dual-hormone 
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
DHparams(1,1) = Fc01; DHparams(2,1) = VdG; DHparams(3,1) = k12;
DHparams(4,1) = Ag; DHparams(5,1) = tmaxG; DHparams(6,1) = EGP0;
%% Insulin Kinetic from Hovorka
tmaxI = 55; Ke = 0.138; VdI = 0.12;
DHparams(7,1) = tmaxI; DHparams(8,1) = Ke; DHparams(9,1) = VdI;
%% Insulin Dynamic
ka1 = 0.006; ka2 = 0.06; ka3 = 0.03;
Sf1 = 51.2e-4; Sf2 = 8.2e-4; Sf3 = 520e-4;
DHparams(10,1) = ka1; DHparams(11,1) = ka2; DHparams(12,1) = ka3; 
DHparams(13,1) = Sf1; DHparams(14,1) = Sf2; DHparams(15,1) = Sf3;
%% Glucagon Kinetic Model (from Lv)
kg1 = 0.0065; kg2 = 0.02777; kge1 = 0.0772; kge2 = 0.0357; VdGG = 0.19;  
DHparams(16,1) = kg1; DHparams(17,1) = kge1;
DHparams(18,1) = kge2; DHparams(19,1) = VdGG; 
DHparams(23,1) = kg2;
%% Glucagon Dynamic Model (Joseph)
SfGG = 0.0171; kc = 0.06; kg3 = 140;
DHparams(20,1) = SfGG; DHparams(21,1) = kc; DHparams(22,1) = kg3;
%% Include the average weight as a model parameter
DHparams(24,1)=weight;

ModPar=DHparams;
end