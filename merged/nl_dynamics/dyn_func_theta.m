function [state_dot] = dyn_func_theta(state, control)

%% Parameters
mass = 1.3474;
Iyy = 0.0667;
g = 9.806;

%%
Va = state(1);
alpha = state(2);
q = state(3);
pitch = state(4);

dE = control(1);


%% Dynamics model
[LIFT, DRAG, M] = lon_aero(Va, alpha, q, dE);

thrust = 0;
b_thrust_ang = 0;
M_thrust = 0;

Va_dot = -DRAG / mass + cos(alpha + b_thrust_ang) * thrust / mass - g * sin(pitch - alpha);
alpha_dot = (-LIFT / mass - sin(alpha + b_thrust_ang) * thrust / mass + g * cos(pitch - alpha) + q * Va) / Va;
q_dot = (M + M_thrust) / Iyy;
pitch_dot = q;

state_dot = [Va_dot; alpha_dot; q_dot; pitch_dot];

end

