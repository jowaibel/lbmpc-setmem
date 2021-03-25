
% Load models
lon_LTI_models;
addpath('nl_dynamics')

% Select model (linear and nonlinear)
mdl = mdls.gamma;
A = mdl.sys.A;
B = mdl.sys.B;
x_trim = mdl.x_trim; nx = length(x_trim);
u_trim = mdl.u_trim; nu = length(u_trim);

dyn_func = @dyn_func_gam;

% Create sim for model
dt = 0.01; % Simulation discretization
sim = SimulatorClass(dyn_func, mdl, dt);

%% Setup Set Membership Estimation
% State constraints; Hx*x <= bx
% Hx=[1 0 0 0;...
%     -1 0 0 0;...
%     0 1 0 0;...
%     0 -1 0 0;...
%     0 0 1 0;...
%     0 0 -1 0;...
%     0 0 0 1;...
%     0 0 0 -1];
% bx=[inf; inf; inf; inf; inf; inf; inf; inf]; %irrelevant at the moment
% 
% % Input constraints; Hu*u <= bu
% Hu=[1;-1];
% bu=[inf;inf];   %irrelevant at the moment
% 
% % Instantiate (Discrete) System
% Ad = eye(size(A,1))+dt*A; Bd = dt*B;
% sysd = LinearSystem(Ad, Bd, Hx, bx, Hu, bu);
% 
% A0 = A; B0 = B(:,1);
% A0(1,1:2) = 1.1 * A0(1,1:2); B0(1) = 1.0 * B0(1);
% 
% % Discrete system
% Ad0 = eye(4) + dt*A0; Bd0 = dt*B0;
% 
% AB = [Ad0, Bd0];
% 
% ABi = zeros(4,5,2);
% ABi(1,1,1) = dt;        % Delta X_u
% ABi(1,2,2) = dt;        % Delta X_w
% 
% % initial bounds on parameters
% uncert = 0.3;
% H_theta = [eye(2); -eye(2)];
% h_theta = uncert * abs([A0(1,1), A0(1,2)])';
% h_theta(3:4) = h_theta(1:2);
% 
% % define initial paramter set
% Omega{1} = Polyhedron(H_theta,h_theta);
% 
% w_max = 0.05;
% % Define additive polytopic uncertainty description
% Hw = [1 0; -1 0; ...
%     0 1; 0 -1];
% hw = w_max * ones(4,1);
% W = Polyhedron(Hw, hw);
% 
% % instantiate set membership estimator
% sm = SetMembership(Omega{1},W, ABi, AB);

%% Simulate
t0 = 0;
t_end = 2;
T = (t0:dt:t_end);

nSteps = ceil(t_end/dt);
T = (t0:dt:t_end);
X = [sim.x_trim, zeros(4, nSteps)]; % nonlinear state trajectory
Xl = X;                     % linear state trajectory

% Possible identification inputs
ident_dt = 0.16; % dt in 3*dt, 2*dt, 1*dt, 1*dt
ident2211 = kron([1 1 -1 -1 1 -1], [ones(1, ident_dt/dt)]);
ident3211 = kron([1 1 1 -1 -1 1 -1], [ones(1, ident_dt/dt)]);
identSig = ident3211;

% Create input trajectory containing identification sequence
% U = zeros(2,length(T));
U = repmat(sim.u_trim, 1, length(T));
idxStart = find(T>0.25, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * identSig;
% idxStart = find(T>2.25, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * identSig;

theta_hat = []; % estimated values

sim.simulate(t0, t_end, mdl.x_trim, U);
sim.plotStateTrajectrory();

for iStep = 1:nSteps
%     X(:,iStep+1)  = sim.simulate_one_step(X(:,iStep), U(:,iStep));
%     [Xl(:,iStep+1), ~] = sim.simulate_one_step_lin(Xl(:,iStep), U(:,iStep));
    
    % Update Set Membership Estimation
%     [Omega{iStep+1}, setD{iStep}] = sm.update(Xl(:,iStep+1), Xl(:,iStep), U(1,iStep));    % update set membership
%     theta_hat = [theta_hat, sm.theta_hat];  % estimate parameter (center of the estimated set)
end

% Plot trajectory
%sim.plotStateTrajectrory(T, X, Xl);

return
%%
figure;
plot(1:1:nSteps, A0(1,1)*ones(1,nSteps)+theta_hat(1,:), 1:1:nSteps, A(1,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("X_u")
figure;
plot(1:1:nSteps, A0(1,2)*ones(1,nSteps)+theta_hat(2,:), 1:1:nSteps, A(1,2)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("X_w")
figure;
plot(1:1:nSteps, A0(1,3)*ones(1,nSteps)+theta_hat(3,:), 1:1:nSteps, A(1,3)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("X_q")
figure;
plot(1:1:nSteps, B0(1,1)*ones(1,nSteps)+theta_hat(4,:), 1:1:nSteps, B(1,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("X_{delta_e}")
figure;
plot(1:1:nSteps, A0(2,1)*ones(1,nSteps)+theta_hat(5,:), 1:1:nSteps, A(2,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("Z_u")
figure;
plot(1:1:nSteps, A0(2,2)*ones(1,nSteps)+theta_hat(6,:), 1:1:nSteps, A(2,2)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("Z_w")
figure;
plot(1:1:nSteps, A0(2,3)*ones(1,nSteps)+theta_hat(7,:), 1:1:nSteps, A(2,3)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("Z_q")
figure;
plot(1:1:nSteps, B0(2,1)*ones(1,nSteps)+theta_hat(8,:), 1:1:nSteps, B(2,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("Z_{delta_e}")
figure;
plot(1:1:nSteps, A0(3,1)*ones(1,nSteps)+theta_hat(9,:), 1:1:nSteps, A(3,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("M_u")
figure;
plot(1:1:nSteps, A0(3,2)*ones(1,nSteps)+theta_hat(10,:), 1:1:nSteps, A(3,2)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("M_w")
figure;
plot(1:1:nSteps, A0(3,3)*ones(1,nSteps)+theta_hat(11,:), 1:1:nSteps, A(3,3)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("M_q")
figure;
plot(1:1:nSteps, B0(3,1)*ones(1,nSteps)+theta_hat(12,:), 1:1:nSteps, B(3,1)*ones(nSteps), '--');
xlabel("timesteps")
ylabel("M_{delta_e}")