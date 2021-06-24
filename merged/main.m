close all
clear
clc
rng('default');

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

% sim.plotStateTrajectory(X');

data.real.data = table(T, U, X, XP);
data.real.dt = dt;
data.real.nx = size(X, 2);
data.real.nu = size(U, 2);
data.real.nSamples = length(data.real.data.X);
n={'V_a'  '\alpha' 'q'  '\Theta'};
un={'m/s' 'rad' 'rad/s' 'rad'};
for i=1:4
    subplot(3,2,i)
    plot(data.real.data.T-data.real.data.T(1),data.real.data.XP(:,i));
    title(n(i))
    ylabel(un(i))
end
subplot(3,1,3)
plot(data.real.data.T-data.real.data.T(1),data.real.data.U',"k")
title("\delta_e");
ylabel("rad");
xlabel("Time (s)");

clear T U X XP dt resampledStateSeqs resampledControlSeqs seqToUse process_noise_abs
clear dyn_func
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
% JA = logical([ % x = [V(a) a(lpha) q theta].
%     1 0 0 0;   % X_V      X_a       0      -g
%     0 1 1 0;   % Z_V/V0   Z_a/V0    1      0
%     0 0 1 0;   % M_V      M_a       M_q    0
%     0 0 0 0]); % 0        0         1      0
% JB = logical([ % u = [elev  (throttle)]
%     0;         % X_el
%     1;         % Z_el/V0
%     1;         % M_el
%     0]);       % 0

J = [JA JB];
smConfig.np = sum(J(:));
smConfig.idxJ = find(J); clear J
% Initial guess for dynamics: XFLR model - but only uncertain values chosen in JA and JB
model0 = mdls.pitch;

% For simulation-generated data, artificially worsen the initial guess
model0.x_trim = mdls.xflr_pitch.x_trim;
model0.u_trim = mdls.xflr_pitch.u_trim;
model0.sys.A(JA(:)) = mdls.xflr_pitch.sys.A(JA(:)); clear JA
model0.sys.B(JB(:)) = mdls.xflr_pitch.sys.B(JB(:)); clear JB
% Select data for estimation
identData = data.real;

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


plotHandles.f1 = figure; smData.iterations = 0; smData.w_bounds = [];
clear dTheta_final_bounds_last

%%
smConfig.recursive_estimation = false;
smConfig.trim_optim = true; 
smConfig.term_crit = 1e-3; % The estimation tries to tighten the dTheta uncertainty bounds until the certainty range in all parameters decreases less than term_crit.
scale_w=1; % term to scale w noise bounds

% Define additive polytopic uncertainty description
smData.DW = identData.data.XP' - smData.AB_ms * [identData.data.X - model0.x_trim', identData.data.U - model0.u_trim(1)']' - model0.x_trim;
smData.w_bound = scale_w*max(abs(smData.DW), [], 2);
% smData.w_highest = 1.5 * smData.w_bound;
% smData.w_lowest = 0.5 * smData.w_bound;

if exist('dTheta_final_bounds_last', 'var'), fprintf('Warmstarting'); end

num=0;
while (true)
    % initialize additive noise
    num=num+1;
    if num==1
        if smConfig.trim_optim
            % initialize trim point
            f=@(trim) sum(sum(identData.data.XP' - trim(1:4,1) - smData.AB_ms * [identData.data.X - trim(1:4,1)', identData.data.U - trim(5,1)']'))^2;
            x0 =  [model0.x_trim; model0.u_trim(1)];
            options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
            [x, fval] = fminunc(f,x0,options);
            model0.x_trim=x(1:sim.nx,1);
            model0.u_trim(1)=x(sim.nx+1,1);
            % compute new w_bounds for new trim estimation
            smData.DW = identData.data.XP' - smData.AB_ms * [identData.data.X - model0.x_trim', identData.data.U - model0.u_trim(1)']' - model0.x_trim;
            smData.w_bound = scale_w*max(abs(smData.DW), [], 2);
        end
           
        % initialize tuples to save data of each iteration
        smData.A_ms = smData.AB_ms(:,1:sim.nx); smData.B_ms = smData.AB_ms(:, sim.nx+1 : sim.nx+sim.nu);
        X_ms = sim.simulate_linear_discr(smData.A_ms, smData.B_ms, model0.x_trim, model0.u_trim(1), identData.data.X(1,:)', identData.data.U(:,1)');   % simulate initial guess and estimated model
        NRMSE_ms_it = [sim.calculateNRMSE(identData.data.X', X_ms(:,1:end-1))];
        trims_it= [model0.x_trim; model0.u_trim(1)];
        AB_it= [smData.AB_ms];
        smData.w_bounds = [smData.w_bound];
        dTheta_bounds=[smData.dTheta_bounds];
        dTheta_hat=[smData.dTheta_hat];
    end
    
    
    % Define initial parameter set
    smData.Omega{1} = Polyhedron(smData.H_theta, smData.h_theta);
    Hh_theta = smData.H_theta .* smData.h_theta;
    smData.dTheta_bounds(1,:) = Hh_theta(logical(repmat(eye(smConfig.np),2,1))); clear Hh_theta
    
    % Define disturbance bounds
    Hw = [eye(sim.nx); -eye(sim.nx)];
    hw = repmat(smData.w_bound, 2, 1);
    smData.W = Polyhedron(Hw, hw); clear Hw hw
    
    % Instantiate set membership estimator
    sm = SetMembership(smData.Omega{1}, smData.W, smData.ABi, smData.AB0);
    
    % Keep track of w_max evolution
    if size(smData.w_bounds,2) - 1 > smData.iterations(end), smData.iterations = [smData.iterations smData.iterations(end)+1]; end
    set(0, 'currentfigure', plotHandles.f1);
    subplot(length(smData.w_bound), 1, 1); plot(smData.iterations, smData.w_bounds(1,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on, title('w\_{bound}')
    subplot(length(smData.w_bound), 1, 2); plot(smData.iterations, smData.w_bounds(2,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on
    subplot(length(smData.w_bound), 1, 3); plot(smData.iterations, smData.w_bounds(3,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on
    subplot(length(smData.w_bound), 1, 4); plot(smData.iterations, smData.w_bounds(4,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(smData.iterations), grid on, drawnow
    
    
    fprintf('\nRunning SM estimation...'); tic;
    if smConfig.recursive_estimation
        disp("to do")
    else
        % Estimate from all samples in one step
        smData.Omega{smConfig.nSteps+1} = sm.update(identData.data.XP', identData.data.X', identData.data.U(:,1)',model0.x_trim, model0.u_trim(1));
        if smData.Omega{smConfig.nSteps+1}.isEmptySet
            % Restart estimation with shrinked W - repeat last step with
            % larger w
            smData.w_bound=1.1*smData.w_bounds(:,end-1);
%             AB_ms = AB_it(:,end-(sim.nx+sim.nu):end)
%             smData.AB_ms= AB_ms;
%             smData.A_ms = AB_ms(:,1:sim.nx); smData.B_ms = AB_ms(:, sim.nx+1 : sim.nx+sim.nu); clear AB_ms
%             model0.x_trim=trims_it(1:sim.nx,end);
%             model0.u_trim(1)=trims_it(sim.nx+1,end);
            fprintf('\nEmpty set, breaking...');
            break
        else
            smData.dTheta_hat(smConfig.nSteps+1, :) = sm.theta_hat';         % estimate parameter (center of the estimated set)
            smData.dTheta_hat(:,:) = interp1([1 2],[smData.dTheta_hat(1,:); smData.dTheta_hat(end,:)], linspace(1, 2, smConfig.nSteps+1));
            
            smData.dTheta_bounds(smConfig.nSteps+1, :) = sm.theta_bounds';
            smData.dTheta_bounds(:,:) = interp1([1 2],[smData.dTheta_bounds(1,:); smData.dTheta_bounds(end,:)], linspace(1, 2, smConfig.nSteps+1));
            
            iStep = smConfig.nSteps; % Set steps to meet following 'estimation completed' condition
            
            
            smData.AB_ms = sm.get_AB();
            
            smData.DW = identData.data.XP' - model0.x_trim - smData.AB_ms * [identData.data.X - model0.x_trim', identData.data.U - model0.u_trim(1)']';
            smData.w_bound = scale_w*max(abs(smData.DW), [], 2);
            
            if smConfig.trim_optim
                f=@(trim) sum(sum(identData.data.XP' - trim(1:sim.nx,1) - smData.AB_ms * [identData.data.X - trim(1:sim.nx,1)', identData.data.U - trim(sim.nx+1,1)']'))^2;
                x0 =  [model0.x_trim; model0.u_trim(1)];
                options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
                [x, fval] = fminunc(f,x0,options);
                model0.x_trim=x(1:sim.nx,1);
                model0.u_trim(1)=x(sim.nx+1,1);
            end
            
            % SAVE iteration results
            trims_it = [trims_it, [model0.x_trim; model0.u_trim(1)]];
            AB_it=[AB_it , smData.AB_ms];
            smData.A_ms = smData.AB_ms(:,1:sim.nx); smData.B_ms = smData.AB_ms(:, sim.nx+1 : sim.nx+sim.nu);
            X_ms = sim.simulate_linear_discr(smData.A_ms, smData.B_ms, model0.x_trim, model0.u_trim(1), identData.data.X(1,:)', identData.data.U(:,1)');   % simulate initial guess and estimated model
            NRMSE_ms_it = [NRMSE_ms_it, sim.calculateNRMSE(identData.data.X', X_ms(:,1:end-1))];
            smData.w_bounds = [smData.w_bounds, smData.w_bound];
            dTheta_bounds=[dTheta_bounds , smData.dTheta_bounds];
            dTheta_hat=[dTheta_hat , smData.dTheta_hat];
            fprintf('\nTighten disturbance set  W.');
        end
    end
    if iStep == smConfig.nSteps % Estimation completed, i.e., all samples used (either recursively or at once)
        if all(abs(smData.w_bounds(:,end)-smData.w_bounds(:,end-1))<smConfig.term_crit) || ... 
                sum(smData.w_bounds(:,end)-smData.w_bounds(:,end-1))>0
            fprintf(['\n---- Estimation complete (' num2str(round(toc,1)) 's). w_bounds span < termination criterion ' ...
                num2str(smConfig.term_crit) '.']);
            clear id du dx dxp iStep
            break
        else
            fprintf(['\n-- All samples used (' num2str(round(toc,1)) 's).' ...
            '----------\n']); %something happend and doesn't work with fprintf
        end
    end
end
clear dTheta_final_bounds_size dTheta_final_bounds_last_size dTheta_final_bounds_size_diff

[v,index]=min(NRMSE_ms_it);
% Set-Membership estimated system - charge ms data at iteration with minimal
% NRMSE value
smData.dTheta_hat=dTheta_hat(:,smConfig.np*(index-1)+1:smConfig.np*index);
smData.dTheta_bounds=dTheta_bounds(:,2*smConfig.np*(index-1)+1:2*smConfig.np*index);
AB_ms = AB_it(:,(sim.nx+sim.nu)*(index-1)+1:(sim.nx+sim.nu)*index)
smData.AB_ms=AB_ms;
smData.A_ms = AB_ms(:,1:sim.nx); smData.B_ms = AB_ms(:, sim.nx+1 : sim.nx+sim.nu); clear AB_ms
model0.x_trim=trims_it(1:sim.nx,index);
model0.u_trim(1)=trims_it(sim.nx+1,index);
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
end, clear iP

disp('Indices changed in A and B:')
disp([smData.A_ms ~= smData.A0, smData.B_ms ~= smData.B0])

% Simulate membership-identified model and plot (data, init model, membership model) trajectories
X_ms = sim.simulate_linear_discr(smData.A_ms, smData.B_ms, model0.x_trim, model0.u_trim(1), identData.data.X(1,:)', identData.data.U(:,1)');   % simulate initial guess and estimated model
sim.plotStateTrajectoryComparison(...
    identData.data.T, identData.data.X', 'data', ...
    identData.data.T, X_init(:,1:end-1), 'init', ...
    identData.data.T, X_ms(:,1:end-1), 'ms');
NRMSE_ms = sim.calculateNRMSE(identData.data.X', X_ms(:,1:end-1));
NRMSE_init = sim.calculateNRMSE(identData.data.X', X_init(:,1:end-1));

fprintf(['\nNormalized RMSE improved from ' num2str(NRMSE_init) ' (initial guess) to ' num2str(NRMSE_ms) ' (estimated model)\n']);
%% 

figure
plot(NRMSE_ms_it)
title("NRMSE set membership over iterations");xlabel('it');ylabel('NRMSE')

pred_1step_ms=model0.x_trim+smData.AB_ms*[identData.data.X'-model0.x_trim;identData.data.U(:,1)'-model0.u_trim(1)];
pred_1step_ini=trims_it(1:sim.nx,1)+AB_it(:,1:(sim.nx+sim.nu))*[identData.data.X'-trims_it(1:sim.nx,1); ...
    identData.data.U(:,1)'-trims_it(sim.nx+1,1)];
figure
sgtitle("Estimations one-step ahead")
for var = 1:sim.nx
    subplot(2,2,var)
    plot(pred_1step_ms(var,:)) %state against next state
    title(char(mdls.xflr_pitch.sys.StateName(var)))
    hold on
    plot(pred_1step_ini(var,:)) %state against next state
    hold on
    plot(identData.data.XP(:,var)')
    legend(["M.S estimation 1-step","ini","data"])
end

% return
%% Maneuver with linearized model and membership model the non linear simulator
see_progress=true; %show progress bar
%AB_ms=[A_ms B_ms];
x_trim_ms = model0.x_trim; 
u_trim_ms = model0.u_trim(1);
x_trim_true = mdls.pitch.x_trim;
u_trim_true = mdls.pitch.u_trim(1);
dt = identData.dt;
DXP=identData.data.XP';
DX=identData.data.X';
DU=identData.data.U';
nx=sim.nx;
nu=sim.nu;
A_ms=smData.A_ms;
B_ms=smData.B_ms;
AB_ms=smData.AB_ms;
A_true = sim.Ade;
B_true = sim.Bde;
AB_true = [A_true, B_true];

% Set up reference
ref_q = [repmat(deg2rad(0),1,40)];
ref_amplitude = deg2rad(20);
ref_freq = 1/4;  % in Hz
ref_endtime = 8;    % in s
ref_theta = repmat(ref_amplitude,1,100);     % Step from trim value to new theta
ref_theta = x_trim_true(4)-ref_amplitude + ref_amplitude* square(2*pi*ref_freq*(0:dt:ref_endtime));      % Square wave reference
ref_theta = x_trim_true(4)-ref_amplitude + ref_amplitude* cos(2*pi*ref_freq*(0:dt:ref_endtime));      % Sine wave reference
ref = [ref_theta];
s = [4]; %state(s) to control in concordance with ref order (in this case applying only theta control)
H = 5; % MPC horizon
opt_steps = size(ref,2)-H;
disp("Controlling state variable ["+char(mdls.xflr_pitch.sys.StateName(s))+"] for "+dt*size(ref,2)+" seconds, horizon H="+H+" steps ("+dt*H+" seconds)")

% Calculate noise bounds:
DW_ms = DXP - x_trim_ms - AB_ms * [DX-x_trim_ms; DU-u_trim_ms];
% DW_true = DXP_true - x_trim_true - AB_true * [DX_true-x_trim_true; DU_true-u_trim_true];
w_max_ms = max(abs(DW_ms), [], 2);
% w_max_true = max(abs(DW_true), [], 2);

%w_max_ms = zeros(4,1);

Hw = [eye(nx);-eye(nx)];
W_ms = Polyhedron(Hw,[w_max_ms; w_max_ms]);
% W_true = Polyhedron(Hw,[w_max_true; w_max_true]);

% Ancillary/Tube Controller for Constraint Tightening
use_tube_controller = true;
if use_tube_controller
    [K_ms, P_ms] = dlqr(smData.A_ms,smData.B_ms,eye(nx),10*eye(nu));     % Q = 10*eye(nx), R = 10*eye(nu)
    K_ms = -K_ms;
    K_true = 0;     % true model has no uncertainty so we don't need constraint tightening
else
    K_ms = 0;
    K_true = 0;
end

% Calculate disturbance reachable sets
upper_ms=[0;0;0;0];
lower_ms=[0;0;0;0];
F_ms{1} = Polyhedron(zeros(0,nx),[]);   %start with empty polyhedron
F_true{1} = Polyhedron(zeros(0,nx),[]);   %start with empty polyhedron
if see_progress progressbar= waitbar(0, 'Starting'); end % Init single bar
for i=1:H
    if use_tube_controller
        F_ms{i+1} = (A_ms+B_ms*K_ms)*F_ms{i}+ W_ms;
        F_true{i+1} = F_true{i};    % all sets empty because true model has no uncertainty
    else
        F_ms{i+1} = plus(A_ms*F_ms{i},W_ms,'vrep');
        F_true{i+e} = F_true{i};    % all sets empty because true model has no uncertainty
    end
    F_ms{i+1}.minHRep;
    F_true{i+1}.minHRep;
    F_ms{i+1}.minVRep;
    F_true{i+1}.minVRep;
    if i>1
        proj_ms=[projection(F_ms{i},[1]).H ; projection(F_ms{i},[2]).H ; projection(F_ms{i},[3]).H ; projection(F_ms{i},[4]).H];
        upper_ms=[upper_ms, abs(proj_ms(1:2:end,2)./proj_ms(1:2:end,1))];
        lower_ms=[lower_ms, -abs(proj_ms(2:2:end,2)./proj_ms(2:2:end,1))];
     
    end
    if see_progress
        waitbar(i/H, progressbar, sprintf('Dist. reachable set Progress: %d %%', floor(i/H*100))); % Update progress bar
    end
end
if see_progress close(progressbar); end

figure
for k =1:nx
    subplot(4,1,k); fill([1:H, fliplr(1:H)],[lower_ms(k,:), fliplr(upper_ms(k,:))],'g');
    if k==1 title("Disturb. reach. sets M.S"); end
    ylabel(char(mdls.xflr_pitch.sys.StateName(k)));
end

%%
% specify states domain
H_x=[eye(nx);-eye(nx)];
% h_x=[13.5 ; 0 ; .4 ; .2; 
%      -12 ; 3 ; .4 ; .2];
h_x = [40; deg2rad(15); 8; deg2rad(20);       
       -5; deg2rad(15); 8; deg2rad(70)]; 

X_set=Polyhedron(H_x,h_x);
for i=1:H
    X_t_ms{i+1} = X_set - F_ms{i+1};
    X_t_true{i+1} = X_set;
    X_t_ms{i+1}=X_t_ms{i+1}.minHRep;
    X_t_true{i+1}=X_t_true{i+1}.minHRep;
end


% Wind disturbance
q_wind = -2;     % in m/s
wind_start = 3; % in seconds
wind_end = 3.1;   % in seconds

%with membership model
x0 = x_trim_true;
X_ms = x0;
U_MPC_ms = [];
u_old = 0;
if see_progress progressbar= waitbar(0, 'Starting'); end
for t = 1:opt_steps
    
    if((t>wind_start/dt)&&(t<wind_end/dt))
        x0(3) = x0(3) + q_wind;
    end
    
    [U,X,u] = MPC_tight(A_ms,B_ms,x0,ref(:,t:end),s,H,F_ms,X_t_ms,K_ms, u_old,x_trim_ms,u_trim_ms);
%     x0=A_true*x0+B_true*u;
%     x0 = sim.simulate_one_step(x0, u);
    x0 = sim.simulate_nonlinear_step(sim.dyn_func, x0, u, dt);
    X_ms=cat(2,X_ms,x0);
    U_MPC_ms=cat(2,U_MPC_ms,u);
    if isnan(u)
        fprintf('\nMPC with est. model infeasible');
        break
    end
    if see_progress
        waitbar(t/opt_steps, progressbar, sprintf('M.S model MPC Progress: %d %%', floor(t/opt_steps*100))); % Update progress bar
    end
    u_old = u;
end
if see_progress close(progressbar); end

%with linearized model
x0 = x_trim_true;
X_true = x0;
U_MPC_true = [];
u_old = 0;
if see_progress progressbar= waitbar(0, 'Starting'); end % Init single bar
for t = 1:opt_steps
    
    if((t>wind_start/dt)&&(t<wind_end/dt))
        x0(3) = x0(3) + q_wind;
    end
    
    [U,X,u]=MPC_tight(A_true,B_true,x0,ref(:,t:end),s,H,F_true,X_t_true,K_true, u_old,x_trim_true,u_trim_true);
%     x0=A_true*x0+B_true*u;
%     x0 = sim.simulate_one_step(x0, u);
    x0 = sim.simulate_nonlinear_step(sim.dyn_func, x0, u, dt);
    X_true=cat(2,X_true,x0);
    U_MPC_true=cat(2,U_MPC_true,u);
    if isnan(u)
        fprintf('\nMPC with true model infeasible');
        break
    end
    if see_progress
        waitbar(t/opt_steps, progressbar, sprintf('True model MPC Progress: %d %%', floor(t/opt_steps*100))); % Update progress bar
    end
    u_old = u;
end
if see_progress close(progressbar); end

%%
time_ms = 0:dt:(size(X_ms,2)-1)*dt;
time_true = 0:dt:(size(X_true,2)-1)*dt;
time = 0:dt:opt_steps*dt;
figure
for k=1:nx
    subplot(4,2,2*k-1); plot(time_ms, X_ms(k,:)); hold on; ylabel([mdls.xflr_pitch.sys.StateName{k}, ' in ', mdls.xflr_pitch.sys.StateUnit{k}]);
    if k==1 title("MPC with estimated model"); end
    plot(time,h_x(k)*ones(opt_steps+1,1),"--"); hold on;
    plot(time,-h_x(k+4)*ones(opt_steps+1,1),"--");
    xlabel('time in s');
end
for k=1:nx
    subplot(4,2,2*k); plot(time_true, X_true(k,:)); hold on; ylabel([mdls.xflr_pitch.sys.StateName{k}, ' in ', mdls.xflr_pitch.sys.StateUnit{k}]);
    if k==1 title("MPC with true model"); end
    plot(time,h_x(k)*ones(opt_steps+1,1),"--"); hold on; 
    plot(time,-h_x(k+4)*ones(opt_steps+1,1),"--")
    xlabel('time in s');
end

figure
for it_s = 1:size(s,1)
    subplot(2+size(s,1),1,it_s)
    plot(time_true, X_true(s(it_s),1:end),"r")
    hold on
    plot(time_ms, X_ms(s(it_s),1:end),"b")
    hold on
    plot(time, ref(it_s, 1:opt_steps+1),"k")
    legend(["True","Estimated","Reference"])
    ylabel([mdls.xflr_pitch.sys.StateName{s(it_s)}, ' in ', mdls.xflr_pitch.sys.StateUnit{s(it_s)}])
    xlabel('time in s');
    title('Reference tracking of controlled state')
end

obj_true=sum((X_true(s,:)-ref(:,1:size(X_true,2))).^2,1);
obj_ms=sum((X_ms(s,:)-ref(:,1:size(X_ms,2))).^2,1);
subplot(2+size(s,1),1,size(s,1)+1)
plot(time, obj_true,"r")
hold on
plot(time, obj_ms,"b")
legend(["True","Estimated"])
title("Objective function value")
xlabel('time in s');

subplot(2+size(s,1),1,size(s,1)+2)
cs_true=cumsum(obj_true);
cs_ms=cumsum(obj_ms);
plot(time, cs_true,"r")
hold on
plot(time, cs_ms,"b")
legend(["True: "+cs_true(end),"Estimated: "+cs_ms(end)])
title("Cumulative objective function value")
xlabel('time in s');

figure;
plot(time(1:end-1), U_MPC_ms, '-', time(1:end-1), U_MPC_true, '--');
legend('Estimated', 'True');
title('Elevator deflection');
xlabel('time in s');
ylabel('eta in rad');