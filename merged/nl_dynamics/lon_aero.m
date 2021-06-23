function [LIFT, DRAG, M] = lon_aero(Va, alpha, q, dE)

b = 1.8;
c = 0.18523;
AR = 10.016;
S = 0.32347;

e_oswald = 0.9;
CD0 = 0.035221;

CL0 = 0.85305;
CLa = 5.6602;
Cm0 = -0.097313;
Cma = -1.1554;

CLq = 0;
Cmq = -30.2134;

CLde = 0;
Cmde = -1.2639;

%%
V0 = 12;
rho = 1.15;

CL = CL0 + CLa * alpha + CLq * c / (2.0 * V0) * q + CLde * dE;
CD = CD0 + CL * CL / (pi * e_oswald * AR);
Cm = Cm0 + Cma * alpha + c / (2.0 * V0) * Cmq * q + Cmde * dE;

dyn_press = 0.5 * rho * Va * Va;
LIFT = dyn_press * S * CL;
DRAG = dyn_press * S * CD;
M = dyn_press * S * c * Cm;

end

