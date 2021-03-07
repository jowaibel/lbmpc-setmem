% Nomenclature
% 
% States:
% vx, (vy), vz  components of velocity v in body frame, also called [u (v) w]
% Va            apparent airspeed
% alpha         angle of attack (AoA), angle of airplane nose above air velocity vector
% q             pitch rate (also omega_y, written as wy), angular velocity about the pitch axis, q = theta_dot
% theta         pitch angle, angle of airplane nose above horizon
% gamma         flight path angle, angle of flight path above horizon
% 
% Note that the vector v = [vx vy vz] is called 'velocity' (oriented) and 
% (capital) V is called 'speed' (unoriented), with V = norm(v).
% 
% Va = norm(v) = norm([vx vz])
% alpha = atan(vz/vx)
% theta = gamma + alpha
% 
% All three state representations represent the same dynamics in
% principle. Some are more suitable for some situations, e.g., the  
% measurements on the real aircraft are vx, vz, q, theta.
% However, Va, alpha are more intuitive to understand the aerodynamics.
% gamma tells you about the efficiency of the flight (how steep do glide 
% descent without motor).
% You can also see that numerically, the trim points are computed slightly
% different. The trim (or steady) state/control is the solution of the condition
% x_dot = 0.
% 
% Control inputs:
% dE            elevator deflection angle
% dF            throttle setting, dimensionless (0 ... 1)
% ! All the model are with dF = 0. Because the thrust model 
% (thrust = f(dF, Va)) is not well identified ! So we could actually just
% use the first column of the B matrix (only one input: dE).

%% x = [Va, alpha, q, theta]
state_names = {'Va'  'alpha' 'q'  'theta'};
state_units = {'m/s' 'rad' 'rad/s' 'rad'};
input_names = {'dE', 'dF'};
input_units = {'rad', '-'};
output_names = state_names;

x_trim = [12.5582, -0.0720904, -4.02425e-09, -0.166115]';
u_trim = [-0.0111063, 0]';

A = [ 
  -0.147491     5.86024           0    -9.76269;
  -0.123808    -9.96186           1   0.0733105;
0.000229934    -94.8472    -19.1422           0;
          0           0           1           0];
B = [ 
           0    -0.425882;
           0 -0.000669201;
    -103.754     0.112586;
           0            0];
      
%% x = [Va, alpha, q, gamma]
state_names = {'Va'  'alpha' 'q'  'gamma'};
state_units = {'m/s' 'rad' 'rad/s' 'rad'};
input_names = {'dE', 'dF'};
input_units = {'rad', '-'};
output_names = state_names;

x_trim = [12.5591, -0.0721008, -1.31188e-08, -0.0940332]';
u_trim = [-0.0110968, 0]';

A = [ 
  -0.147494    -3.90244           0    -9.76268;
  -0.123791     -9.8892           1   0.0733121;
0.000229919    -94.8597    -19.1447           0;
   0.123791      9.8892           0  -0.0733121];
B = [ 
           0    -0.425916;
           0 -0.000669563;
    -103.768     0.112595;
           0  0.000669563];

%% x = [vx, vz, q, theta]
state_names = {'vx'  'vz' 'q'  'theta'};
state_units = {'m/s' 'm/s' 'rad/s' 'rad'};
input_names = {'dE', 'dF'};
input_units = {'rad', '-'};
output_names = state_names;

x_trim = [12.4573, -0.888532, 1.38132e-08, -0.164542]';
u_trim = [-0.0119153, 0]';

A = [ 
-0.274288 -0.217785  0.888532  -9.67355;
  -2.2507   -9.7801   12.4573   1.60623;
-0.534138   -7.4919  -18.9315         0;
        0         0         1         0];
B = [ 
         0  -0.422556;
         0 -0.0221452;
  -102.612   0.111838;
         0          0];


%% Create a state space model and induce step
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

sys_ss = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
canon(sys_ss)

% transfer function
tfsys = tf(sys_ss);

figure
% step response (dE --> all state variables)
[Y, T] = step(tfsys);
subplot(4,1,1)
plot(T, Y(:,1))
title([sys_ss.StateName(1)])
grid on

subplot(4,1,2)
plot(T, Y(:,2))
title([sys_ss.StateName(2)])
grid on

subplot(4,1,3)
plot(T, Y(:,3))
title([sys_ss.StateName(3)])
grid on

subplot(4,1,4)
plot(T, Y(:,4))
title([sys_ss.StateName(4)])
grid on