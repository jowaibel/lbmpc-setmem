state_names = {'vx', 'vy', 'vz', 'wx', 'wy', 'wz', 'roll', 'pitch'};
input_names = {'dF', 'dE', 'dR', 'dA'};
output_names = {'vx', 'vy', 'vz', 'wx', 'wy', 'wz', 'roll', 'pitch'};

% % % ------------------------------------------------------------------------------------
% % % Initial model (from CFD)
% % x_trim = [13.32, 0, -0.57, 0, 0, 0, -0.0781, 0];
% % 
% % A = [...
% %      -0.116107          0  0.0398224          0   0.468065          0          0   -9.77611;
% %              0  -0.596423          0  -0.786757          0   -12.9347    9.77611          0;
% %       -1.90471          0   -9.89957         -0    11.8868          0          0    0.76507;
% %              0   -5.11791         -0   -26.4363         -0    5.73717          0          0;
% %      -0.307941          0   -7.72014         -0   -8.86108         -0          0          0;
% %              0    3.21013         -0   -2.48935         -0   -2.23854          0          0;
% %              0          0          0          1         -0 -0.0782592          0          0;
% %              0          0          0          0          1         -0         -0          0];
% %          
% %          
% % B = [... 
% %     -0.48332        0        0;
% %            0  5.00952        0;
% %     -6.79563        0        0;
% %            0  7.30501 -228.495;
% %     -83.6201        0        0;
% %            0 -38.0764 -8.40074;
% %            0        0        0;
% %            0        0        0];
% ------------------------------------------------------------------------------------
% Best identified model (closest to 'true')
x_trim = [13.8, 0, -1.19, 0, 0, 0, 0, -0.1948, 0];
u_trim = [0, 0, 0]; % (Throttle = 0)

A = [... 
-0.318435         0 -0.496135         0      1.19         0         0  -9.62053
        0 -0.967949         0     -1.19         0  -11.9268   9.62053         0;
 -2.33071         0  -10.8187        -0      13.8         0         0   1.89815;
        0  -4.57177        -0  -18.4337        -0   9.28813         0         0;
-0.685885         0  -8.30197        -0  -20.1746        -0         0         0;
        0   2.70627        -0  -1.49769        -0  -2.23426         0         0;
        0         0         0         1        -0 -0.197302         0         0;
        0         0         0         0         1        -0        -0         0];
B = [... 
-0.478507         0         0         0;
        0         0    10.644         0;
        0         0         0         0;
        0         0    16.825   -129.41;
        0  -126.219         0         0;
        0         0  -32.0859  -2.37201;
        0         0         0         0;
        0         0         0         0];
% ------------------------------------------------------------------------------------

%% Reduce for longitudinal dynamics
x_idx = [1 3 5 8];
u_idx = [1 2];
state_names = state_names(x_idx);
input_names = input_names(u_idx);
output_names = output_names(x_idx);
A = A(x_idx, x_idx);
B = B(x_idx, u_idx);


%% =======================================================================
%%% Enter Longitudinal Model directly

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

%% x = [vx, vz, q, gamma]
state_names = {'vx'  'vz' 'q'  'pitch'};
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