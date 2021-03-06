clc
clear

% Load models
mdls = lon_LTI_models();
addpath('nl_dynamics')

% Select model (linear and nonlinear)
mdl = mdls.pitch_withPos; % Select linear model

% Remove throttle input
mdl.sys = ss(mdl.sys.A, mdl.sys.B(:,1), mdl.sys.C, mdl.sys.D(:,1), 'StateName', mdl.sys.StateName, 'StateUnit', mdl.sys.StateUnit, ...
    'InputName', mdl.sys.InputName{1}, 'InputUnit', mdl.sys.InputUnit{1}, ...
     'OutputName', mdl.sys.OutputName);
mdl.u_trim = mdl.u_trim(1);
    
Ac = mdl.sys.A
Bc = mdl.sys.B
x_trim = mdl.x_trim; nx = length(x_trim);
u_trim = mdl.u_trim; nu = length(u_trim);

dyn_func = @dyn_func_theta_withPos; % Select (same) nonlinear model

% Create sim for model
dt = 0.01; % Simulation discretization
sim = SimulatorClass(dyn_func, mdl, dt);

% Simulate
t0    = 0;
t_end = 4;
T = (t0:dt:t_end);       % Time vector (s)

% Possible identification inputs
ident_dt = 0.16; % dt in 3*dt, 2*dt, 1*dt, 1*dt
ident2211 = kron([1 1 -1 -1 1 -1], [ones(1, ident_dt/dt)]);
ident3211 = kron([1 1 1 -1 -1 1 -1], [ones(1, ident_dt/dt)]);
identSig = ident3211;

% Create input trajectory containing identification sequence
% U = zeros(2,length(T));
U = repmat(u_trim, 1, length(T));
idxStart = find(T>0.1, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * identSig; % Insert ident command starting at 0.25s
idxStart = find(T>2.25, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * -identSig;

% Simulate both nonlinear and linear model in parallel (same initial state and control trajectory)
% x0 = mdl.x_trim;
% sim.simulate(t0, t_end, x0, U);
% sim.plotStateTrajectrory();

%% Simulate
t0 = 0;
t_end = 2;
T = (t0:dt:t_end);

nSteps = ceil(t_end/dt);
T = (t0:dt:t_end);
X = [sim.x_trim, zeros(6, nSteps)]; % nonlinear state trajectory
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

%theta_hat = []; % estimated values

x0 = mdl.x_trim;
sim.simulate(t0, t_end, x0, U);
sim.plotStateTrajectrory();
return
for iStep = 1:nSteps
    X(:,iStep+1)  = sim.simulate_one_step(X(:,iStep), U(:,iStep));
%     [Xl(:,iStep+1), ~] = sim.simulate_one_step_lin(Xl(:,iStep), U(:,iStep));
    
    % Update Set Membership Estimation
%     [Omega{iStep+1}, setD{iStep}] = sm.update(Xl(:,iStep+1), Xl(:,iStep), U(1,iStep));    % update set membership
%     theta_hat = [theta_hat, sm.theta_hat];  % estimate parameter (center of the estimated set)
end

% Plot trajectory
sim.plotStateTrajectrory(T, X);

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