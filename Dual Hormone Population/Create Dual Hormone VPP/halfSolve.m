  function F = halfSolve(u,Gsp,HovPar)

% first, for the 'guessed' u (insulin infusion), solve for steady-state
%   state variables

  x0 = [ones(8,1)*0.1;zeros(5,1)];
  for t = 1:7000
      dxdt = hovssJEY(x0,u,HovPar);
      x0 = x0 + dxdt;
      x0(x0 < 0) = 0;
  end
  Vg = HovPar(2);
  GmM = x0(1)/Vg;
  F = Gsp - GmM;