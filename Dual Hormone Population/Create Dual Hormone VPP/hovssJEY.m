function F = hovssJEY(x,u,HovPar)

mass_patient = HovPar(24); % mass in kg
F01 = HovPar(1);      % nominal value, F01 = 0.0097 mmol / kg-min
Vg  = HovPar(2);      % nominal value, Vg = 0.16 L / kg
k12 = HovPar(3);      % nominal value, k12 = 0.066 min^-1
Ag = HovPar(4);      
tmaxG = HovPar(5);    
EGP0 = HovPar(6);     % nominal value, EGP0 = 0.0161 mmol / kg-min
tmaxI = HovPar(7);    % nominal value, k = 0.67 dimensionless (MKK had 0.57)
Ke = HovPar(8);       % nominal value, ka1i = 0.0112 min^-1
VdI = HovPar(9);     % nominal value, ka2i = 0.0210 min^-1
ka1 = HovPar(10);     % nominal value, ka1 = 0.006 min^-1
ka2 = HovPar(11);     % nominal value, ka2 = 0.06 min^-1
ka3 = HovPar(12);     % nominal value, ka3 = 0.03 min^-1
SfIT = HovPar(13);    % nominal value, SfIT = 51.2e-4 (mU/L)^-1.(min^-1)
SfID = HovPar(14);    % nominal value, SfID = 8.2e-4 (mU/L)^-1.(min^-1)
SfIE = HovPar(15);    % nominal value, SfIE = 520e-4 (mU/L)^-1
kg1 = HovPar(16);     % kg1 = 0.011124 min^-1
kge1 = HovPar(17);    % kge1 = 0.13227 min^-1
kge2 = HovPar(18);    % kge1 = 0.14608 min^-1
Vgg = HovPar(19);     % Vgg = 0.0385 L/kg
Sfgg = HovPar(20);    % Sfgg = 0.0052 (ng/L)^-1
kc = HovPar(21);      % kc = 0.01 min^-1
kg3 = HovPar(22); 
kg2 = HovPar(23);     % kg2 = 0.02777 min^-1

%
%  evaluate calculated parameters
kb1 = SfIT*ka1;   % (mU/L)^-1.(min^-2)
kb2 = SfID*ka2;   % (mU/L)^-1.(min^-2)
kb3 = SfIE*ka3;   % (mU/L)^-1.(min^-1)
kg = Sfgg*kc; % (ng/L)^-1.(min^-1)
%
Ug = 0; % no meal
%
Q1  = x(1);             % mmol/kg
Q2  = x(2);             % mmol/kg
S1 = x(3);             
S2 = x(4);             
I = x(5);             
x1  = x(6);             % min^-1
x2  = x(7);             % min^-1
x3  = x(8);             % unitless
Q1g = x(9);             % mg/kg
Q2g = x(10);            % mg/kg
Q3g = x(11);            % mg/kg
Y = x(12);              
Z = x(13);              

%
G = Q1/Vg;            % (mmol/kg)/(L/kg) == mmol/L
                      % for mg/dL, multiply mmol/L by 18
Gg = Q3g * 1E6/Vgg;   % (mg/kg)/(L/kg) == mg/L * 1000000 = ng/L(or pg/ml)

if G >= 4.5
    F01_c = F01;
else
    F01_c = F01 * G/4.5;
end
%
if G >= 9
    Fr = 0.003 * (G-9) * Vg;
else
    Fr = 0;
end
% constrain EGP -- added 4/1/2014
EGP = EGP0*(1 - x3 + Y + kg3*Z);
if EGP < 0
    EGP = 0;
end

uInsulin = u;  % mU/kg/min
uGlucagon = 0; % mg/kg/min - No glucagon during steady-state for insulin
% Glucose compartment derivatives
dQ1dt  = -((F01_c/(Vg*G)) + x1)*Q1 + k12*Q2 - Fr + Ug + EGP;
dQ2dt  = x1*Q1 - (k12 + x2)*Q2;
% Insulin subcutaneous compartment derivatives
dS1dt = uInsulin - S1/tmaxI;  % u = insulin infusion, mU/kg/min
dS2dt = S1/tmaxI - S2/tmaxI;
dIdt = S2/(tmaxI*VdI) - Ke*I;
% Insulin action compartment derivatives
dx1dt  = -ka1*x1 + kb1*I;
dx2dt  = -ka2*x2 + kb2*I;
dx3dt  = -ka3*x3 + kb3*I;
% Glucagon subcutaneous compartment derivatives
dQ1gdt = uGlucagon - kg1*Q1g - kge1*Q1g;
dQ2gdt = kg1*Q1g - kg2*Q2g;
dQ3gdt = kg2*Q2g - kge2*Q3g;
% Glucagon action compartment derivatives
dYdt = -kc*Y + kg*Gg;
dZdt = kg*kg2*Q2g - kg*kge2*Q3g - kc*Z;
F = [dQ1dt;dQ2dt;dS1dt;dS2dt;dIdt;dx1dt;dx2dt;dx3dt;dQ1gdt;dQ2gdt;dQ3gdt;dYdt;dZdt];

