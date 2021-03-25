clear
clc
     
Ac_uw = [ 
-0.274288 -0.217785  0.888532  -9.67355;
  -2.2507   -9.7801   12.4573   1.60623;
-0.534138   -7.4919  -18.9315         0;
        0         0         1         0];
Bc_uw = [ 
         0;
         0;
  -102.612;
         0];
% Discretized system
Ts = 0.01;   % =time_res;
Ad_uw = eye(4) + Ts*Ac_uw;
Bd_uw = Ts*Bc_uw;    

% State constraints; Hx*x <= bx
Hx=[1 0 0 0;...
    -1 0 0 0;...
    0 1 0 0;...
    0 -1 0 0;...
    0 0 1 0;...
    0 0 -1 0;...
    0 0 0 1;...
    0 0 0 -1];
bx=[inf; inf; inf; inf; inf; inf; inf; inf]; %irrelevant at the moment

% Input constraints; Hu*u <= bu
Hu=[1;-1];
bu=[inf;inf];   %irrelevant at the moment

% Instantiate System
sys_uw = LinearSystem(Ac_uw, Bc_uw, Hx, bx, Hu, bu);
clear A_nom_uw B_nom_uw Hx bx Hu bu
     
%% Simulate

Ts = 0.01;
t_end = 10;
traj.t = 0:Ts:t_end';
traj.x = zeros(4, length(traj.t));
traj.u = zeros(1, length(traj.t));

ident_dt = 0.16;
ident2211 = kron([1 1 -1 -1 1 -1], [ones(1, ident_dt/Ts)]);
traj.u(:, 50:50-1+length(ident2211)) = deg2rad(5) * ident2211;

x_trim = [12.4573, -0.888532, 1.38132e-08, -0.164542]';
u_trim = [-0.0119153, 0]';

theta_hat = []; % estimated values

for i=1:length(traj.t)-1
    traj.x(:,i+1) = sys_uw.A * (traj.x(:,i)) + sys_uw.B * (traj.u(:,i));   % simulate one step
    
    %[Omega{i+1},D{i}] = sm.update(traj.x(:,i+1), traj.x(:,i), traj.u(:,i));    % update set membership
    
    %theta_hat = [theta_hat, sm.theta_hat];  % estimate parameter (center of the estimated set)
end
traj.x = traj.x + repmat(x_trim, 1, length(traj.t));

figure;
plot(traj.t, traj.x(4,:))
