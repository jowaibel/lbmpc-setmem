close all
clear
clc

% Load linear and nonlinear model
mdls = lon_LTI_models();
addpath('nl_dynamics')

% ## USER ## Select model (linear and nonlinear)
mdl = mdls.pitch; % Select linear model

% Remove throttle input
mdl.sys = ss(mdl.sys.A, mdl.sys.B(:,1), mdl.sys.C, mdl.sys.D(:,1), ...
    'StateName', mdl.sys.StateName, 'StateUnit', mdl.sys.StateUnit, ...
    'InputName', mdl.sys.InputName{1}, 'InputUnit', mdl.sys.InputUnit{1}, ...
    'OutputName', mdl.sys.OutputName);
mdl.u_trim = mdl.u_trim(1);

dyn_func = @dyn_func_theta; % ## USER ## Select (same) nonlinear model

%% Generate simulation data / ground truth = best identified nonlinear model

% Create sim for model
dt = 0.01; % Simulation discretization
process_noise_abs = [1 0.5 3 0] * dt; % Maximum value of absolute process noise
sim = SimulatorClass(dyn_func, mdl, dt, process_noise_abs);

clear dt process_noise_abs

% Simulate
t0    = 0;
t_end = 4;
T = (t0 : sim.dt : t_end+sim.dt);       % Time vector (s)
clear t0 t_end

% Possible identification inputs
ident_dt = 0.16; % dt in 3*dt, 2*dt, 1*dt, 1*dt
ident2211 = kron([1 1 -1 -1 1 -1], [ones(1, ident_dt/sim.dt)]);
ident3211 = kron([1 1 1 -1 -1 1 -1], [ones(1, ident_dt/sim.dt)]);
identSig = ident3211;
identSig = identSig .* (rand(1,length(identSig)) * 0.2 + 0.9); % Randomize a bit to gain more information
clear ident_dt ident2211 ident3211

% Create input trajectory containing identification sequence
% U = zeros(2,length(T));
U = repmat(sim.u_trim, 1, length(T));
idxStart = find(T>0.1, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * identSig; % Insert ident command starting at 0.25s
idxStart = find(T>2, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * -identSig;
clear idxStart identSig

% Simulate both nonlinear and linear model in parallel (same initial state and control trajectory)
sim.simulate_all(T, mdl.x_trim, U);

% Plot trajectories with different nonlin/lin / discretization settings
% Plot: {nonlinear, linearized} continuous dynamics, integrated with ode45
%       discretized linearized {d: Euler, de: exact}
sim.plotStateTrajectory();

% Assign variables to comply with existing framework
T = T(1:end-1)';
U = U(1:end-1)';

X = sim.nonlin_state_traj(:,1:end-1)';
XP = sim.nonlin_state_traj(:,2:end)';
data.nonlin.data = table(T, U, X, XP);
data.nonlin.dt = T(2) - T(1);
data.nonlin.nx = size(X, 2);
data.nonlin.nu = size(U, 2);

X = sim.linc_state_traj(:,1:end-1)';
XP = sim.linc_state_traj(:,2:end)';
data.lin.data = table(T, U, X, XP);
data.lin.dt = T(2) - T(1);
data.lin.nx = size(X, 2);
data.lin.nu = size(U, 2);
clear T U X XP

% Cut whole time and input signals
cutIdx = (data.nonlin.data.T >= 0); %& data.nonlin.data.T <= 2);
data.nonlin.data = data.nonlin.data(cutIdx,:);
data.nonlin.nSamples = length(data.nonlin.data.X);

data.lin.data = data.lin.data(cutIdx,:);
data.lin.nSamples = length(data.lin.data.X);

clear cutIdx

%% Load real data
load('ident_data') % real data
seqToUse = 2; % 2-6 are suitable
U = resampledControlSeqs{seqToUse}(:, {'elevator'}).Variables;
X = resampledStateSeqs{seqToUse}(:, {'airspeed', 'alpha1', 'pitchRate_smooth', 'pitch'}).Variables;
XP = X;
T = seconds(resampledStateSeqs{seqToUse}.Time);

% Shrink dataset by one to obtain XP
U = -U(1:end-1, :);
X = X(1:end-1, :);
XP = XP(2:end, :);
T = T(1:end-1);

dt = T(2) - T(1);

% Simulate with best nonlinear model and plot mismatch with data
sim = SimulatorClass(dyn_func, mdl, dt);
sim.simulate_all(T, X(1,:)', U(:,1)');
%sim.plotStateTrajectory(X');

data.real.data = table(T, U, X, XP);
data.real.dt = dt;
data.real.nx = size(X, 2);
data.real.nu = size(U, 2);
data.real.nSamples = length(data.real.data.X);

clear T U X XP dt resampledStateSeqs resampledControlSeqs seqToUse process_noise_abs
clear dyn_func

% data.real = data.real(cutIdx,:);
% nSamples_real = size(data.real.X, 1);
% ------------------------------------------------------------------------

%% Configure set membership estimation

% Initial guess for dynamics: XFLR model
JA = logical([ % x = [V(a) a(lpha) q theta].
    1 1 0 0;   % X_V      X_a       0      -g
    1 1 0 0;   % Z_V/V0   Z_a/V0    1      0
    1 1 1 0;   % M_V      M_a       M_q    0
    0 0 0 0]); % 0        0         1      0
JB = logical([ % u = [elev  (throttle)]
    1;         % X_el
    1;         % Z_el/V0
    1;         % M_el
    0]);       % 0

JA = logical([ % x = [V(a) a(lpha) q theta].
    0 0 0 0;   % X_V      X_a       0      -g
    0 1 0 0;   % Z_V/V0   Z_a/V0    1      0
    0 1 1 0;   % M_V      M_a       M_q    0
    0 0 0 0]); % 0        0         1      0
JB = logical([ % u = [elev  (throttle)]
    0;         % X_el
    1;         % Z_el/V0
    1;         % M_el
    0]);       % 0

J = [JA JB];
smConfig.np = sum(J(:));
smConfig.idxJ = find(J); clear J

% Initial guess for dynamics: XFLR model - but only uncertain values chosen in JA and JB
model0 = mdls.pitch;

% Select data for estimation
identData = data.lin;

% For simulation-generated data, artificially worsen the initial guess
model0.x_trim = mdls.xflr_pitch.x_trim;
model0.u_trim = mdls.xflr_pitch.u_trim;
model0.sys.A(JA(:)) = mdls.xflr_pitch.sys.A(JA(:)); clear JA
model0.sys.B(JB(:)) = mdls.xflr_pitch.sys.B(JB(:)); clear JB

% Discretize initial model
Ac0 = model0.sys.A;
Bc0 = model0.sys.B(:,1);

% smData.A0 = eye(sim.nx) + identData.dt * Ac0; % Euler approximation
% smData.B0 = identData.dt * Bc0;
M = expm([Ac0 Bc0; zeros(sim.nu, sim.nx+sim.nu)]*identData.dt); % Exact discretization with matrix-exponential-function
smData.A0 = M(1:sim.nx, 1:sim.nx);
smData.B0 = M(1:sim.nx, (sim.nx+1):(sim.nx+sim.nu)); clear M

smData.AB0 = [smData.A0 smData.B0]; clear Ac0 Bc0 A0 B0
disp('AB0 = '), disp(smData.AB0)

% Simulate and plot initial model over ground truth
X_init = sim.simulate_linear_discr(smData.A0, smData.B0, model0.x_trim, model0.u_trim(1), identData.data.X(1,:)', identData.data.U(:,1)');

% Prepare set membership
smData.ABi = zeros([size(smData.AB0) smConfig.np]);
for iP = 1:smConfig.np
    ABi_ = zeros(size(smData.AB0));
    ABi_(smConfig.idxJ(iP)) = identData.dt;
    smData.ABi(:,:,iP) = ABi_;
end, clear ABi_ iP

% Bounds on delta parameters
smData.H_theta = [eye(smConfig.np); -eye(smConfig.np)];
theta_rel_uncert = 50.0;
% Ac0Bc0 = [Ac0, Bc0];
% Ac0Bc0_pos = Ac0Bc0; Ac0Bc0_pos(Ac0Bc0 < 0) = 0;                % ### Asymmetric uncertainty
% Ac0Bc0_neg = Ac0Bc0; Ac0Bc0_neg(Ac0Bc0 > 0) = 0;                % ### Asymmetric uncertainty
% h_theta = theta_uncert * [Ac0Bc0_pos(idxJ); -Ac0Bc0_neg(idxJ)]; % ### Asymmetric uncertainty
% % AB0_pos = AB0; AB0_pos(AB0 < 0) = 0;                % ### Asymmetric uncertainty (Discrete system)
% % AB0_neg = AB0; AB0_neg(AB0 > 0) = 0;                % ### Asymmetric uncertainty
% % h_theta = ( [-AB0_neg(idxJ); AB0_pos(idxJ)] + theta_rel_uncert * [AB0_pos(idxJ); -AB0_neg(idxJ)] ) ./ dt; % ### Bounds to avoid sign flip + uncertainty range in same sign
smData.h_theta = repmat(theta_rel_uncert * abs( smData.AB0(smConfig.idxJ) ) ./ identData.dt, 2, 1); % ### Symmetric uncertainty
clear theta_rel_uncert

% Set number of ident samples
smConfig.nSteps = identData.nSamples;
% smConfig.nSteps = 160; % Manual override

% Setup sample drawing method from data
%identConfig.dataIdx = randperm(smData.nSteps);                                                                   % Random selection from dataset
%identConfig.dataIdx = 1:floor(smData.nSteps/smData.nSamples):smData.nSteps; dataIdx = dataIdx(1:smData.nSteps);  % Equally distributed selection over time (= downsampling)
smConfig.dataIdx = 1:1:smConfig.nSteps;                                                                         % first nSteps points from dataset

smConfig.nSteps = length(smConfig.dataIdx);
smConfig.stepIdx = 0:smConfig.nSteps;

smData.AB_ms = smData.AB0;
smData.dTheta_hat = zeros(smConfig.nSteps+1, smConfig.np);       % estimated values
smData.dTheta_bounds = zeros(smConfig.nSteps+1, 2*smConfig.np);
smData.setD = cell(1, smConfig.nSteps+1);

% Define additive polytopic uncertainty description
smData.DW = identData.data.XP' - smData.AB_ms * [identData.data.X - model0.x_trim', identData.data.U - model0.u_trim(1)']' - model0.x_trim;
smData.w_bound = max(abs(smData.DW), [], 2);
smData.w_highest = 1.5 * smData.w_bound;
smData.w_lowest = 0.5 * smData.w_bound;

plotHandles.f1 = figure; smData.iterations = 0; smData.w_bounds = [];
clear dTheta_final_bounds_last

%%
smConfig.recursive_estimation = true;

smConfig.term_crit = 5; % The estimation tries to tighten the dTheta uncertainty bounds until the certainty range in all parameters decreases less than term_crit.

if exist('dTheta_final_bounds_last', 'var'), fprintf('Warmstarting'); end

while (true)
    % Define initial parameter set
    smData.Omega{1} = Polyhedron(smData.H_theta, smData.h_theta);
    Hh_theta = smData.H_theta .* smData.h_theta;
    smData.dTheta_bounds(1,:) = Hh_theta(logical(repmat(eye(smConfig.np),2,1))); clear Hh_theta
    
    % Define disturbance bounds
    Hw = [eye(sim.nx); -eye(sim.nx)];
    smData.w_bound = (smData.w_lowest + smData.w_highest)/2;
    hw = repmat(smData.w_bound, 2, 1);
    smData.W = Polyhedron(Hw, hw); clear Hw hw
    
    % Instantiate set membership estimator
    sm = SetMembership(smData.Omega{1}, smData.W, smData.ABi, smData.AB0);
    
    % Keep track of w_max evolution
    smData.w_bounds = [smData.w_bounds smData.w_bound];
    if size(smData.w_bounds,2) - 1 > smData.iterations(end), smData.iterations = [smData.iterations smData.iterations(end)+1]; end
    set(0, 'currentfigure', plotHandles.f1);
    subplot(length(smData.w_bound), 1, 1); plot(smData.iterations, smData.w_bounds(1,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on, title('w\_{max}')
    subplot(length(smData.w_bound), 1, 2); plot(smData.iterations, smData.w_bounds(2,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on
    subplot(length(smData.w_bound), 1, 3); plot(smData.iterations, smData.w_bounds(3,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on
    subplot(length(smData.w_bound), 1, 4); plot(smData.iterations, smData.w_bounds(4,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on, drawnow
    
    
    fprintf('\nRunning SM estimation...'); tic;
    if smConfig.recursive_estimation
        % Estimate recursively, sample by sample
        iPlot = 1;
        for iStep = smConfig.stepIdx(2:end)
            iSample = smConfig.dataIdx(iStep);
            
            % update set membership
            [smData.Omega{iStep+1}, smData.setD{iStep+1}] = sm.update(identData.data.XP(iSample,:)', identData.data.X(iSample,:)', identData.data.U(iSample,:)');
            
            if smData.Omega{iStep+1}.isEmptySet
                % Restart estimation with larger W
                smData.w_lowest = max(smData.w_lowest, smData.w_bound);
                fprintf([' -- Step ' num2str(iStep) '/' num2str(smConfig.nSteps) ': Empty set, enlarge disturbance set W.']);
                break;
            end
            
            smData.dTheta_hat(iStep+1,:) = sm.theta_hat';  % estimate parameter (center of the estimated set)
            smData.dTheta_bounds(iStep+1,:) = sm.theta_bounds';
            
            if (smConfig.np == 2) %&& (iStep >= nSteps - 4) % Plot 2D set (if estimating only two parameters)
                if ~smData.setD{iStep+1}.isBounded
                    fprintf([' - D not bounded.']);
                    smData.setD{iStep+1} = Polyhedron('A', smData.setD{iStep+1}.A, 'b', smData.setD{iStep+1}.b, 'lb', -1000*ones(2,1), 'ub', 1000*ones(2,1));
                end
                subplot(2,2,iPlot)
                hold off
                plot(smData.setD{iStep+1},'color','b','alpha',0.1)   % Unfalsified set from this measurement
                hold on
                plot(smData.Omega{iStep+1},'alpha',0.1)              % Present param set
                plot(smData.Omega{iStep+1})                        % Intersection
                title(['k=', int2str(iStep)])
                iPlot = iPlot + 1;
                if iPlot > 4, iPlot = 1; end
            end
        end
    else
        % Estimate from all samples in one step
        smData.Omega{smConfig.nSteps+1} = sm.update(identData.data.XP', identData.data.X', identData.data.U');
        if smData.Omega{smConfig.nSteps+1}.isEmptySet
            % Restart estimation with shrinked W
            smData.w_lowest = max(smData.w_lowest, smData.w_bound);
            fprintf('\nEmpty set, enlarge disturbance set W.');
            continue
        else
            smData.dTheta_hat(smConfig.nSteps+1, :) = sm.theta_hat';         % estimate parameter (center of the estimated set)
            smData.dTheta_hat(:,:) = interp1([1 2],[smData.dTheta_hat(1,:); smData.dTheta_hat(end,:)], linspace(1, 2, smConfig.nSteps+1));
            
            smData.dTheta_bounds(smConfig.nSteps+1, :) = sm.theta_bounds';
            smData.dTheta_bounds(:,:) = interp1([1 2],[smData.dTheta_bounds(1,:); smData.dTheta_bounds(end,:)], linspace(1, 2, smConfig.nSteps+1));
            
            iStep = smConfig.nSteps; % Set steps to meet following 'estimation completed' condition
        end
    end
    
    if iStep == smConfig.nSteps % Estimation completed, i.e., all samples used (either recursively or at once)
        smData.AB_ms = sm.get_AB(); % Save current AB estimate
        % Prepare restart of estimation with updated W (warmstart)
        % Restart estimation with shrinked W
        smData.w_highest = min(smData.w_highest, smData.w_bound);
        
        if exist('dTheta_final_bounds_last', 'var')
            % Previous result is available for comparison
            
            % Range between dTheta bounds
            dTheta_final_bounds_size = abs(diff(reshape(smData.dTheta_bounds(end,:), 2, smConfig.np)));
            dTheta_final_bounds_last_size = abs(diff(reshape(dTheta_final_bounds_last, 2, smConfig.np)));
            
            % Change in that range w.r.t. to previous complete estimation
            dTheta_final_bounds_size_diff = dTheta_final_bounds_last_size - dTheta_final_bounds_size;
            if all(abs(dTheta_final_bounds_size_diff) >= -0.1) && ... % There was an improvement &&
                    all(dTheta_final_bounds_size_diff < smConfig.term_crit) % Improvement is small enough to stop
                fprintf(['\nEstimation complete (' num2str(round(toc,1)) 's). dTheta_bounds span < termination criterion ' num2str(smConfig.term_crit) '.' ...
                    '\nFor more accurate estimation, run section again and/or decrease term_crit.\n ==========\n']);
                clear id du dx dxp iStep
                dTheta_final_bounds_last = smData.dTheta_bounds(end,:);
                break
            end
        end
        fprintf([' -- All samples used (' num2str(round(toc,1)) 's). Tighten disturbance set W.\n ----------']);
        dTheta_final_bounds_last = smData.dTheta_bounds(end,:);
    end
end
clear dTheta_final_bounds_size dTheta_final_bounds_last_size dTheta_final_bounds_size_diff

% Set-Membership estimated system
AB_ms = smData.AB_ms
smData.A_ms = AB_ms(:,1:sim.nx); smData.B_ms = AB_ms(:, sim.nx+1 : sim.nx+sim.nu); clear AB_ms

% Output true system (discrete) for comparison
AB_true = [sim.Ade sim.Bde]

% Plot delta and absolute parameter values
for iP = 1:smConfig.np
    figure(iP+10), %clf;
    subplot(1,2,1) % Delta parameter estimation
    hold on;
    patch([smConfig.stepIdx'; flipud(smConfig.stepIdx')], [smData.dTheta_bounds(:,2*iP-1); flipud(smData.dTheta_bounds(:,2*iP))], 'k', 'FaceAlpha', 0.1); % Prediction intervals
    plot(0:smConfig.nSteps, smData.dTheta_hat(:,iP), 'b');
    ylim([mean(smData.dTheta_bounds(:,2*iP)), mean(smData.dTheta_bounds(:,2*iP-1))])
    title(['Param ' num2str(iP) ' dTheta hat'])
    
    subplot(1,2,2) % Parameter estimation
    hold on
    % True value (linearized dynamics) dashed line
    plot(smConfig.stepIdx, repmat(AB_true(smConfig.idxJ(iP)), smConfig.nSteps+1), 'k--');
    % Estimation history
    patch([smConfig.stepIdx'; flipud(smConfig.stepIdx')], smData.AB0(smConfig.idxJ(iP)) + identData.dt * [smData.dTheta_bounds(:,2*iP-1); flipud(smData.dTheta_bounds(:,2*iP))], 'k', 'FaceAlpha', 0.1); % Prediction intervals
    plot(smConfig.stepIdx, smData.AB0(smConfig.idxJ(iP)) + identData.dt * smData.dTheta_hat(:,iP), 'b');
    ylim(smData.AB0(smConfig.idxJ(iP)) + identData.dt * [mean(smData.dTheta_bounds(:,2*iP)), mean(smData.dTheta_bounds(:,2*iP-1))])
    title(['Param ' num2str(iP) ' value'])
end, clear iP AB_true

disp('Indices changed in A and B:')
disp([smData.A_ms ~= smData.A0, smData.B_ms ~= smData.B0])

% Simulate membership-identified model and plot (data, init model, membership model) trajectories
X_ms = sim.simulate_linear_discr(smData.A_ms, smData.B_ms, model0.x_trim, model0.u_trim(1), identData.data.X(1,:)', identData.data.U(:,1)');   % simulate initial guess and estimated model
sim.plotStateTrajectoryComparison(...
    identData.data.T, identData.data.X', 'data', ...
    identData.data.T, X_init(:,1:end-1), 'init', ...
    identData.data.T, X_ms(:,1:end-1), 'ms');

NRMSE_init = sim.calculateNRMSE(identData.data.X', X_init(:,1:end-1));
NRMSE_ms = sim.calculateNRMSE(identData.data.X', X_ms(:,1:end-1));
fprintf(['\nNormalized RMSE improved from ' num2str(NRMSE_init) ' (initial guess) to ' num2str(NRMSE_ms) ' (estimated model)\n']);
return

% These membership estimations should'nt converge, as the initial guess is the horizontal line, its a different plot that the one used before.



figure
sgtitle("Nextstate v/s state, Linearized model (red), data (blue), M.S model (green)")
states=["Va","alpha","q","gamma"];
for var = 1:4
    subplot(2,2,var)
    plot(df_s(:,var),df_ns(:,var), 'bp', 'MarkerSize', 0.5) %state against next state
    hold on
    plot(df_s(:,var),df_lin_ns(:,var), 'rp', 'MarkerSize', 0.5)
    hold on
    plot(df_s(:,var),df_ms_ns(:,var), 'gp', 'MarkerSize', 1)
    xlim([min(df_s(:,var)),max(df_s(:,var))]);
    ylim([min(df_ns(:,var)),max(df_ns(:,var))]);
    xlabel(states(var));
end

figure
sgtitle(["Nextstate(data) - Nextstate(linearized model) [blue]"; ...
    "and Nextstate(data) - Nextstate(membership model) [green]"])
for var = 1:4
    subplot(2,2,var)
    plot(df_s(:,var),df_ns(:,var)-df_lin_ns(:,var), 'bp', 'MarkerSize', 1)
    hold on;
    plot(df_s(:,var),df_ns(:,var)-df_ms_ns(:,var), 'gp', 'MarkerSize', 1)
    xlim([min(df_s(:,var)),max(df_s(:,var))]);
    xlabel(states(var));
end

%% Maneuver with linearized model and membership model the non linear simulator
% In constrast with the membership model, the linearized model can not be obtained without aircraft parameter knowledgement. There are two nominal MPC functions implemented, as the linearized model works with xdot=A(x-xtrim).... and the membership with xdot=Ax+Bu...
% Define maneuver (this maneuver is a stepwise trajectory tracking maneuver
% looking to align alpha and gamma):
alpha_ref=[repmat(deg2rad(-4),1,200) repmat(deg2rad(0),1,500)];
gamma_ref=alpha_ref;

%with linearized model
x_ini=x_trim';
x0=x_ini';
X_lin=x0;
for t = 1:400
    [U,X,u]=MPC_nom_lin(x0,dt,alpha_ref(:,t:end),gamma_ref(:,t:end),20);
    x0=x0+dyn_func_gam(x0, u)*dt;
    X_lin=cat(2,X_lin,x0);
    if isnan(u)
        break
    end
end

%with membership model
x_ini=x_trim';
x0=x_ini';
X_ms=x0;
for t = 1:400
    [U,X,u]=MPC_nom(A_ms,B_ms,x0,dt,alpha_ref(:,t:end),gamma_ref(:,t:end),20);
    x0=x0+dyn_func_gam(x0, u)*dt;
    X_ms=cat(2,X_ms,x0);
    if isnan(u)
        break
    end
end

figure
subplot(3,1,1)
plot(X_lin(2,:),"r")
hold on
plot(X_ms(2,:),"g")
hold on
plot(alpha_ref,"k")
legend(["linearized","M.S","ref."])
ylabel("alpha")
subplot(3,1,2)
plot(X_lin(4,:),"r")
hold on
plot(X_ms(4,:),"g")
hold on
plot(gamma_ref,"k")
legend(["linearized","M.S","ref."])
ylabel("gamma")

figure
sgtitle("Objective function evaluation")
plot((X_lin(2,:)-alpha_ref(1,1:size(X_lin,2))).^2+(X_lin(4,:)-gamma_ref(1,1:size(X_lin,2))).^2,"r")
hold on
plot((X_ms(2,:)-alpha_ref(1,1:size(X_ms,2))).^2+(X_ms(4,:)-gamma_ref(1,1:size(X_ms,2))).^2,"g")
legend(["linearized","M.S"])

figure
sgtitle("Cumulative objective function evaluation")
plot(cumsum((X_lin(2,:)-alpha_ref(1,1:size(X_lin,2))).^2+(X_lin(4,:)-gamma_ref(1,1:size(X_lin,2))).^2),"r")
hold on
plot(cumsum((X_ms(2,:)-alpha_ref(1,1:size(X_ms,2))).^2+(X_ms(4,:)-gamma_ref(1,1:size(X_ms,2))).^2),"g")
legend(["linearized","M.S"])

