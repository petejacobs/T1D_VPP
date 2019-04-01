function [Ad,Bd,Cd,Dd] = Plant(X,ModPar)
% This function implements the glucoregulatory model for state X and model
% parameters ModPar.  The outputs of the model are in the form
% x(t+1) = Ad x(t) + Bd u(t)  + Dd
% y(t+1) = Cd x(t+1)
% See paper Resalat et al. 2019, PLOS One for further details
% Inputs:
% X:  States of the model
% ModPar: Parameters of the model
% Outputs:
% 
Ts = 5;
Fc01 = ModPar(1);     % nominal value, Fc01 = 0.0097 mmol / kg-min
Vdg  = ModPar(2);     % nominal value, Vdg = 0.16 L / kg
k12 = ModPar(3);      % nominal value, k12 = 0.066 min^-1
Ag = ModPar(4);       % nominal value, Ag = 0.8 
tmaxG = ModPar(5);    % nominal value, tmaxG = 40
EGP0 = ModPar(6);     % nominal value, EGP0 = 0.0161 mmol / kg-min
tmaxI = ModPar(7);     
Ke = ModPar(8);      
VdI = ModPar(9);       
ka1 = ModPar(10);     % nominal value, ka1 = 0.006 min^-1
ka2 = ModPar(11);     % nominal value, ka2 = 0.06 min^-1
ka3 = ModPar(12);     % nominal value, ka3 = 0.03 min^-1
Sf1 = ModPar(13);     % nominal value, SfIT = 51.2e-4
Sf2 = ModPar(14);     % nominal value, SfID = 8.2e-4
Sf3 = ModPar(15);     % nominal value, SfIE = 520e-4
% Sf1 = M_E_PIU * M_E_PGU * Sf1; Sf2 = M_E_PIU * M_E_PGU * Sf2;
% Sf3 = M_E_HGP * Sf3;

kb1 = Sf1*ka1;        % L / mU
kb2 = Sf2*ka2;        % L / mU
kb3 = Sf3*ka3;        % L / mU - min^-1

Q1  = X(1);           % mmol/L
Q2  = X(2);
S1 = X(3);
S2 = X(4);
I = X(5);
x1  = X(6);
x2  = X(7);
x3  = X(8);

Num_States = 8; No_Cmplx = Num_States;
G = Q1/Vdg;            % (mmol/kg)/(L/kg) == mmol/L
                       % for mg/dL, multiply mmol/L by 18
                       
Ap = zeros(Num_States,Num_States);
Bp = zeros(Num_States,1);
Dp = zeros(Num_States,1);
Cp = zeros(1,Num_States);

if G < 4.5
    KcO1 = Fc01/(4.5*Vdg); KR = 0; KD = 0;
elseif G >= 4.5 && G < 9
    KcO1 = 0; KR = 0; KD = Fc01;
else
    KcO1 = 0; KR = 0.003; KD = Fc01 - 0.003*9*Vdg;
end

Ap(1,1) = -x1 - KcO1 - KR;
Ap(1,2) = k12; Ap(1,6) = -Q1; Ap(1,8) = -EGP0;
Ap(2,1) = x1; Ap(2,2) = -(k12+x2); Ap(2,6) = Q1; Ap(2,7) = -Q2;

Ap(3,3) = -1/tmaxI;
Ap(4,3) = 1/tmaxI; Ap(4,4) = -1/tmaxI;
Ap(5,4) = 1/(tmaxI*VdI); Ap(5,5) = -Ke;

Ap(6,6) = -ka1; Ap(7,7) = -ka2; Ap(8,8) = -ka3;
Ap(6,5) = kb1; Ap(7,5) = kb2; Ap(8,5) = kb3;

Bp(3,1) = 1;      

Dp(1,1) = x1*Q1 - KD + EGP0;
Dp(2,1) = -x1*Q1 + x2*Q2;

Cp(1,1)=1;

Delta_t = Ts;
% Dx = zeros(1,1); [Ad, Bd, Cd, ~] = c2dm(Ap,Bp,Cp,Dx,Delta_t);


s = expm([[Ap Bp]*Ts; zeros(1,9)]);
Ad = s(1:No_Cmplx,1:No_Cmplx);
Bd = s(1:No_Cmplx,No_Cmplx+1:end);

Cd = Cp; 
Dd = Dp*Delta_t;


