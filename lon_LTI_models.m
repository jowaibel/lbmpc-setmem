function mdls = lon_LTI_models()
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

% xflr ----------------------------------------------
x_trim = [12.5021, -0.0327415, 8.59739e-08, -0.0680528];
u_trim = [-0.0105831, 0];

A = [ 
 -0.0562547     6.09266  -0.0446645    -9.79989;
  -0.125394    -9.34381    0.887762   0.0276905;
0.000230966    -90.5807     -8.6571           0;
          0           0           1           0];
B = [ 
  -0.190614   -0.423593;
  -0.478997 0.000664789;
    -73.532    0.111981;
          0           0];
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

mdls.xflr_pitch.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.xflr_pitch.x_trim = x_trim;
mdls.xflr_pitch.u_trim = u_trim;
      
% xflr-Pvw-YR ----------------------------------------------
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
       
mdls.pitch.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.pitch.x_trim = x_trim;
mdls.pitch.u_trim = u_trim;

% xflr-Pvw-YR withPos ----------------------------------------------
state_names = {'Va'  'alpha' 'q'  'theta', 'posHor', 'posDown'};
state_units = {'m/s' 'rad' 'rad/s' 'rad', 'm', 'm'};
output_names = state_names;

x_trim = [12.5582, -0.0720904, -4.02425e-09, -0.166115, 0, -100]';
u_trim = [-0.0111063, 0]';

A = [A zeros(4, 2);
    zeros(2, 6)];
B = [B
    0 0;
    0 0];
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

mdls.pitch_withPos.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.pitch_withPos.x_trim = x_trim;
mdls.pitch_withPos.u_trim = u_trim;

%% x = [Va, alpha, q, gamma]
state_names = {'Va'  'alpha' 'q'  'gamma'};
state_units = {'m/s' 'rad' 'rad/s' 'rad'};
input_names = {'dE', 'dF'};
input_units = {'rad', '-'};
output_names = state_names;

% xflr ----------------------------------------------
x_trim = [12.5394, -0.0332745, 7.27206e-09, -0.035336];
u_trim = [-0.0099264, 0];

A = [
-0.0561262   -3.70722 -0.0446645   -9.79988;
  -0.12465   -9.34388   0.887427  0.0276276;
0.00023028   -91.1213   -8.70877          0;
   0.12465    9.34388   0.112573 -0.0276276];
B = [ 
   -0.190614    -0.425118;
   -0.480425  0.000647123;
    -73.9708     0.112383;
    0.480425 -0.000647123];

nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

mdls.xflr_gamma.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.xflr_gamma.x_trim = x_trim;
mdls.xflr_gamma.u_trim = u_trim;

% xflr-Pvw-YR ----------------------------------------------
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

mdls.gamma.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.gamma.x_trim = x_trim;
mdls.gamma.u_trim = u_trim;

% xflr-Pvw-YR with Pos ----------------------------------------------
state_names = {'Va'  'alpha' 'q'  'gamma', 'posHor', 'posDown'};
state_units = {'m/s' 'rad' 'rad/s' 'rad', 'm', 'm'};
output_names = state_names;

x_trim = [12.5591, -0.0721008, -1.31188e-08, -0.0940332, 0, -100]';
A = [A zeros(4, 2);
     zeros(2, 6)];
B = [B
     0 0;
     0 0];
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

mdls.gamma_withPos.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.gamma_withPos.x_trim = x_trim;
mdls.gamma_withPos.u_trim = u_trim;

%% x = [vx, vz, q, theta]
state_names = {'vx'  'vz' 'q'  'theta'};
state_units = {'m/s' 'm/s' 'rad/s' 'rad'};
input_names = {'dE', 'dF'};
input_units = {'rad', '-'};
output_names = state_names;

% xflr ----------------------------------------------
x_trim = [12.4573, -0.888532, 1.38132e-08, -0.164542]';
u_trim = [-0.0119153, 0]';

A = [ 
-0.101401  0.186261  0.317147  -9.78338;
  -1.8702  -9.29274   11.0883  0.665618;
-0.235922  -7.23673  -8.64601         0;
        0         0         1         0];
B = [ 
 -0.385589  -0.422767
  -5.97142 -0.0221563
  -73.4378   0.111894
         0          0];
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);

mdls.xflr_uw.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.xflr_uw.x_trim = x_trim;
mdls.xflr_uw.u_trim = u_trim;

% xflr-Pvw-YR ----------------------------------------------
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
     
mdls.uw.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.uw.x_trim = x_trim;
mdls.uw.u_trim = u_trim;

% xflr-Pvw-YR withPos ----------------------------------------------
state_names = {'vx'  'vz' 'q'  'theta', 'posHor', 'posDown'};
state_units = {'m/s' 'm/s' 'rad/s' 'rad', 'm', 'm'};
output_names = state_names;

x_trim = [12.4573, -0.888532, 1.38132e-08, -0.164542 0 -100]';
u_trim = [-0.0119153, 0]';

A = [A zeros(4, 2);
     zeros(2, 6)];
B = [B
     0 0;
     0 0];
nx = size(A,2);
nu = size(B,2);
C = eye(nx);
D = zeros(nx, nu);
 
mdls.uw_withPos.sys = ss(A, B, C, D, 'StateName', state_names, 'StateUnit', state_units, ...
    'InputName', input_names, 'InputUnit', input_units, ...
    'OutputName', output_names);
mdls.uw_withPos.x_trim = x_trim;
mdls.uw_withPos.u_trim = u_trim;

end
