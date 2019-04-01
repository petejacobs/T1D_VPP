function [Ad,Bd,Cd,Dd] = SH_Glucoregulatory(X, ModPar, Ts, M_E_PIU, M_E_PGU, M_E_HGP)

Fc01 = ModPar(1);
Vdg  = ModPar(2);
k12 = ModPar(3);
Ag = ModPar(4);
tmaxG = ModPar(5);
EGP0 = ModPar(6);
tmaxI = ModPar(7);
Ke = ModPar(8);
VdI = ModPar(9);
ka1 = ModPar(10);
ka2 = ModPar(11);
ka3 = ModPar(12);
Sf1 = ModPar(13);
Sf2 = ModPar(14);
Sf3 = ModPar(15);
Sf1 = M_E_PIU * M_E_PGU * Sf1; Sf2 = M_E_PIU * M_E_PGU * Sf2;
Sf3 = M_E_HGP * Sf3;

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

Num_States = 8;
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
Cd = Cp; Dd = Dp*Delta_t;

%% c2d ZOH method
% s = expm([[Ap Bp]*Ts; zeros(1,Num_States+1)]);
% Ad = s(1:Num_States,1:Num_States);
% Bd = s(1:Num_States,Num_States+1:end);

%% c2d forward method
Ad = Ap*Delta_t + eye(Num_States, Num_States); Bd = Bp*Delta_t; 

