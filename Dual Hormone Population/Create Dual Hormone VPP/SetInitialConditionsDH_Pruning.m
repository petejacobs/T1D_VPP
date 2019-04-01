function [xm, iBas] = SetInitialConditionsDH_Pruning(ModPar)

direc = 0; % 1 means increase, 0 means decrease
iBas = 0.6;
du = iBas/2;
dG = 1;

while abs(dG) > 0.001
    dG = halfSolve(iBas,100/18,ModPar);
    if dG > 0.001
        if direc == 1
            direc = 0;
            du = du/2;
        end
        iBas = iBas - du;
    elseif dG < -0.001
        if direc == 0
            direc = 1;
            du = du/2;
        end
        iBas = iBas + du;
    end
end

x0 = [ones(8,1)*0.1;zeros(5,1)];
for t = 1:7000
    dxdt = hovssJEY(x0,iBas,ModPar);
    x0 = x0 + dxdt;
    x0(x0 < 0) = 0;
    G(t) = x0(1)*18/ModPar(2);
end
xss = x0;


Q1 = xss(1,1);                  % (mmol/kg) - Measurable glucose compartment
Q2 = xss(2,1);                  % (mmol/kg) - Unmeasurable glucose compartment
S1 = xss(3,1);              
S2 = xss(4,1);               
I = xss(5,1);                
x1 = xss(6,1);                 % (min^-1) - Rate constant compartment for glucose distribution (Q1 to Q2)
x2 = xss(7,1);                % (min^-1) - Rate constant compartment for insulin mediated glucose utilization (elimination from Q2)
x3 = xss(8,1);    
x1g = 0;
x2g = 0;
x3g = 0;
Y = 0;
Z = 0;

xm = [Q1;Q2;S1;S2;I;x1;x2;x3;x1g;x2g;x3g;Y;Z];
end