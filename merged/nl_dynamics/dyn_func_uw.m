function [state_dot] = dyn_func_uw(state, control)

%% Parameters
mass = 1.3474;
Iyy = 0.0667;
g = 9.806;

%%
vx = state(1);
vz = state(2);
q = state(3);
theta = state(4);

dE = control(1);


%% Dynamics model
Va = norm([vx, vz]);
alpha = atan(vz / vx);
[LIFT, DRAG, M] = lon_aero(Va, alpha, q, dE);

b_F_thrust = [0 0 0];
M_thrust = 0;

X = -cos(alpha) * DRAG + sin(alpha) * LIFT;
Z = -cos(alpha) * LIFT - sin(alpha) * DRAG;
vx_dot = (X + b_F_thrust(1)) / mass - g * sin(theta) - q * vz;
vz_dot = (Z + b_F_thrust(3)) / mass + g * cos(theta) + q * vx;
q_dot = (M + M_thrust) / Iyy;
theta_dot = q;

state_dot = [vx_dot; vz_dot; q_dot; theta_dot];

end

