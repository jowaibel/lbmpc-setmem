%% Wrapping up everything
% Collect data from simulation
% Initially we have a non-linear model for aircraft flight simulation, from
% which we can generate some data for random states and random actions.
% Let's consider states and actions around the x_trim point to generate
% this data.
close all
clear
clc

% Load models
mdls = lon_LTI_models();
addpath('nl_dynamics')

% ## USER ## Select model (linear and nonlinear)
mdl = mdls.uw; % Select linear model

% Remove throttle input
mdl.sys = ss(mdl.sys.A, mdl.sys.B(:,1), mdl.sys.C, mdl.sys.D(:,1), ...
    'StateName', mdl.sys.StateName, 'StateUnit', mdl.sys.StateUnit, ...
    'InputName', mdl.sys.InputName{1}, 'InputUnit', mdl.sys.InputUnit{1}, ...
    'OutputName', mdl.sys.OutputName);
mdl.u_trim = mdl.u_trim(1);

Ac = mdl.sys.A
Bc = mdl.sys.B
x_trim = mdl.x_trim; nx = length(x_trim);
u_trim = mdl.u_trim; nu = length(u_trim);

dyn_func = @dyn_func_uw; % ## USER ## Select (same) nonlinear model

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
identSig = identSig .* (rand(1,length(identSig)) * 0.2 + 0.9); % Randomize a bit to gain more information

% Create input trajectory containing identification sequence
% U = zeros(2,length(T));
U = repmat(u_trim, 1, length(T));
idxStart = find(T>0.1, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * identSig; % Insert ident command starting at 0.25s
idxStart = find(T>2, 1); U(1, idxStart:idxStart-1+length(identSig)) = deg2rad(20) * -identSig;

% Simulate both nonlinear and linear model in parallel (same initial state and control trajectory)
sim.simulate(t0, t_end, mdl.x_trim, U);

% Plot: {nonlinear, linearized} continuous dynamics, integrated with ode45
%       discretized linearized {d: Euler, de: exact}
sim.plotStateTrajectory();

% Assign variables to comply with existing framework
df_u = U(:,1:end-1)';

df_s = sim.state_traj(:,1:end-1)';
df_ns = sim.state_traj(:,2:end)';

df_lin_s = sim.lin_state_traj(:,1:end-1)';
df_lin_ns = sim.lin_state_traj(:,2:end)';

%obtain linearized approximation of next states (Using A(x-x_trim)+B(u-u_trim))
%df_lin_ns = df_s(:,1:4) + dt * (Ac(:,:)*(df_s(:,1:4)-repmat(x_trim',nSteps,1))'+Bc(:,:)*(df_u(:,1)-repmat(u_trim(1)',nSteps,1))')';


% Check discretization uncertainty
A_true = sim.Ad; % Discretization from truncated Talor series (Euler discretization)
B_true = sim.Bd;
% M = expm([Ac Bc; zeros(nu, nx+nu)]*dt);   % Discretization with matrix-exponential-function
% A_true(:,:) = M(1:nx, 1:nx);
% B_true(:,:) = M(1:nx, (nx+1):(nx+nu));

% Cut whole time and input signals
cutIdx = (T > 0 & T < 1.9);
T = T(cutIdx);
%t0 = T(1);
%t_end = T(end);
%U = U(cutIdx);
df_u = df_u(cutIdx,:);
df_s = df_s(cutIdx,:);
df_ns = df_ns(cutIdx,:);
df_lin_s = df_lin_s(cutIdx,:);
df_lin_ns = df_lin_ns(cutIdx,:);

nData = size(df_s, 1);
df_disc_ns = repmat(x_trim',nData,1) + (A_true*(df_lin_s(:,1:4)-repmat(x_trim',nData,1))' + B_true*(df_u(:,1)-repmat(u_trim(1)',nData,1))')';
disc_uncert = max((df_lin_ns-df_disc_ns),[],1);

%%
% figure
% sgtitle("Nextstate v/s state, Linearized model (red) and data (blue)")
% for var = 1:4
%     subplot(2,2,var)
%     plot(df_s(:,var),df_ns(:,var), 'bp', 'MarkerSize', 1) %state against next state
%     hold on
%     plot(df_s(:,var),df_lin_ns(:,var), 'rp', 'MarkerSize', 1)
%     xlim([min(df_s(:,var)), max(df_s(:,var))]);
%     ylim([min(df_ns(:,var)),max(df_ns(:,var))]);
%     xlabel(mdl.sys.StateName{var});
% end
% figure
% sgtitle("Nextstate(data) - Nextstate(linearized model)")
% for var = 1:4
%     subplot(2,2,var)
%     plot(df_s(:,var),df_ns(:,var)-df_lin_ns(:,var), 'p', 'MarkerSize', 1)
%     xlim([min(df_s(:,var)),max(df_s(:,var))]);
%     ylim([min(df_ns(:,var)-df_lin_ns(:,var)),max(df_ns(:,var)-df_lin_ns(:,var))]);
%     xlabel(mdl.sys.StateName{var});
% end

%% Set Membership Estimation for dimensional derivatives
% I have used the model for xdot=A(x-x_trim)+B*(u-u_trim) as
% xdot=A*x + B*u plus noise to set the initial guess (remember that the
% result from m.s. identification is Ax+Bu). This initial guess can be
% changed with Waibel's XFLR. I have expanded the membership estimation so
% all values in the A and B matrices can be considered uncertain if its
% index is in the J set, so adapt to another models can be made easily.
% This J set to indicate indexes where is uncertainty, has to be though as
% the A and B matrices were merged into one matrix and then flatten, it
% follows the numeration used by Tilman previously in h_theta but
% generalizes for different elements in A and B with uncertainty.

%% Initial guess for dynamics: XFLR model
JA = logical([
    1 1 1 1;
    1 1 1 1;
    1 1 1 0;
    0 0 1 0]);
JB = logical([
    0;
    0;
    0;
    0]);
J = [JA JB];
idxJ = find(J);

% Initial guess for dynamics: XFLR model - but only uncertain values chosen in JA and JB
init_model = mdls.uw;
init_model.sys.A(JA(:)) = mdls.xflr_uw.sys.A(JA(:));
init_model.sys.B(JB(:)) = mdls.xflr_uw.sys.B(JB(:));
%init_model = mdls.gamma;
Ac0 = init_model.sys.A;
Bc0 = init_model.sys.B(:,1);

theta_uncert_true = [abs(Ac0 - Ac)./Ac0, abs(Bc0 - Bc)./Bc0];    % true uncertainty (used for debugging/selecting initial theta_uncert value)
theta_uncert_true(isnan(theta_uncert_true)) = 0

% Discretize initial model
A0 = eye(nx) + dt * Ac0;                      % Euler approximation
B0 = dt * Bc0;
% M = expm([Ac0 Bc0; zeros(nu, nx+nu)]*dt);   % Exact discretization with matrix-exponential-function
% A0(:,:) = M(1:nx, 1:nx);
% B0(:,:) = M(1:nx, (nx+1):(nx+nu));
AB0 = [A0 B0]


% this is to indicate to membership function which parameters has uncertainty
np = sum(J(:));
ABi = zeros([size(AB0) np]);
for iP = 1:np
    ABi_ = zeros(size(AB0));
    ABi_(idxJ(iP)) = dt;
    ABi(:,:,iP) = ABi_;
end, clear ABi_

% Bounds on delta parameters
H_theta=[eye(np); -eye(np)];
theta_uncert = 50.0;
Ac0Bc0 = [Ac0, Bc0];
h_theta = repmat(theta_uncert * abs(Ac0Bc0(idxJ)), 2, 1);


% Get Set Membership estimation
nSteps = nData;
nSteps = 100;
% Select nSteps samples from dataset, either randomly or equally distributed
%dataIdx = randperm(nData);          % Random selection from dataset
%dataIdx = 1:floor(nSteps/nSteps):nSteps; dataIdx = dataIdx(1:nSteps);  % Equally distributed selection over time (= downsampling)
dataIdx = 1:1:nSteps; % first nSteps points from dataset

nSteps = length(dataIdx);
stepIdx = 0:nSteps;

dTheta_hat = zeros(nSteps+1, np); % estimated values
dTheta_bounds = zeros(nSteps+1, 2*np);
setD = cell(1, nSteps+1);

% Define additive polytopic uncertainty description
w_max = 0.1; %0.041;
Hw = [eye(nx); -eye(nx)];

f1 = figure; iterations = 0; w_maxes = [];
clear dTheta_final_bounds_last
%%
recursive_estimation = true;
term_crit = 5; % The estimation tries to tighten the dTheta uncertainty bounds until the certainty range in all parameters decreases less than term_crit.
        
if exist('dTheta_final_bounds_last', 'var'), fprintf('Warmstarting'); end

while (true)
    fprintf(['\nStarting SM estimation with w_max = ' num2str(w_max)]); tic;
    
    w_maxes = [w_maxes w_max];
    if length(w_maxes)-1 > iterations(end), iterations = [iterations iterations(end)+1]; end
    set(0, 'currentfigure', f1);
    plot(iterations, w_maxes, '.-k', 'MarkerSize', 12), title('w\_{max}'), xlabel('it'), xticks(iterations), grid on,  drawnow
    
    % define initial parameter set
    Omega{1} = Polyhedron(H_theta, h_theta);
    dTheta_bounds(1,:) = nonzeros(H_theta .* h_theta);
    
    hw = w_max * ones(2*nx, 1);
    W = Polyhedron(Hw, hw);
    % instantiate set membership estimator
    sm = SetMembership(Omega{1}, W, ABi, AB0);
    
    if recursive_estimation
        % Estimate recursively, sample by sample
        
        % Use nSteps samples for set membership identification
        iPlot = 1;
        for iStep = stepIdx(2:end)
            id = dataIdx(iStep);
            du = df_u(id,1)' - u_trim;
            dx = df_lin_s(id,:)' - x_trim;
            dxp = df_lin_ns(id,:)' - x_trim;
            [Omega{iStep+1}, setD{iStep+1}] = sm.update(dxp, dx, du);    % update set membership
            
            if Omega{iStep+1}.isEmptySet
                % Restart estimation with larger W
                w_max = 1.1 * w_max;
                hw = w_max * ones(2*nx, 1);
                fprintf([' -- Step ' num2str(iStep) '/' num2str(nSteps) ': Empty set, enlarge disturbance set W.']);
                break;
            end
            
            dTheta_hat(iStep+1,:) = sm.theta_hat';  % estimate parameter (center of the estimated set)
            dTheta_bounds(iStep+1,:) = sm.theta_bounds';
            
            if (np == 2) %&& (iStep >= nSteps - 4) % Plot 2D set (if estimating only two parameters)
                if ~setD{iStep+1}.isBounded
                    fprintf([' - D not bounded.']);
                    setD{iStep+1} = Polyhedron('A',setD{iStep+1}.A, 'b', setD{iStep+1}.b, 'lb', -1000*ones(2,1), 'ub', 1000*ones(2,1));
                end
                subplot(2,2,iPlot)
                hold off
                plot(setD{iStep+1},'color','b','alpha',0.1)   % Unfalsified set from this measurement
                hold on
                plot(Omega{iStep+1},'alpha',0.1)              % Present param set
                plot(Omega{iStep+1})                        % Intersection
                title(['k=', int2str(iStep)])
                iPlot = iPlot + 1;
                if iPlot > 4, iPlot = 1; end
            end
        end
    else
        % Estimate from all samples in one step
        DU = df_u(dataIdx,1)' - u_trim;
        DX = df_lin_s(dataIdx,:)' - x_trim;
        DXP = df_lin_ns(dataIdx,:)' - x_trim;
        Omega{nSteps+1} = sm.update_all(DXP, DX, DU);
        if Omega{nSteps+1}.isEmptySet
            % Restart estimation with larger W
            w_max = 1.1 * w_max;
            hw = w_max * ones(2*nx, 1);
            fprintf('\nAll samples used. Empty set, enlarge disturbance set W.');
            break;
        end
        dTheta_hat(nSteps+1,:) = sm.theta_hat';  % estimate parameter (center of the estimated set)
        dTheta_hat(:,:) = interp1([1 2],[dTheta_hat(1,:);dTheta_hat(end,:)], linspace(1,2,nSteps+1));

        dTheta_bounds(nSteps+1,:) = sm.theta_bounds';
        dTheta_bounds(:,:) = interp1([1 2],[dTheta_bounds(1,:);dTheta_bounds(end,:)], linspace(1,2,nSteps+1));
        iStep = nSteps;
    end
    
    if iStep == nSteps
        % Estimation has completed through all samples
        
        % Compute new w_max, in case of having another loop (or warmstart)
        w_max = 0.9 * w_max;
        hw = w_max * ones(2*nx, 1);
        
        if exist('dTheta_final_bounds_last', 'var')
            % Range between dTheta bounds
            dTheta_final_bounds_size = abs(diff(reshape(dTheta_bounds(end,:), 2, 12)));
            dTheta_final_bounds_last_size = abs(diff(reshape(dTheta_final_bounds_last, 2, 12)));
            
            % Change in that range w.r.t. to previous complete estimation
            dTheta_final_bounds_size_diff  = dTheta_final_bounds_last_size - dTheta_final_bounds_size;
            if all(dTheta_final_bounds_size_diff < term_crit)
                % Improvement w.r.t. last cmoplete estimation is
                % sufficiently small, terminate
                fprintf(['\n ==== \nEstimation complete (' num2str(round(toc,1)) 's) with termination criterion ' num2str(term_crit) ', w_max = ' num2str(W.b(1)) '.' ...
                    '\nFor more accurate estimation, run section again and/or decrease term_crit.\n\n']);
                clear id du dx dxp iStep
                dTheta_final_bounds_last = dTheta_bounds(end,:);
                break
            end
        end
        fprintf(['\nAll samples used (' num2str(round(toc,1)) 's). Tighten disturbance set W.\n ----']);
        dTheta_final_bounds_last = dTheta_bounds(end,:);
    end
end
%Omega_end = Omega{nSteps+1}; %figure, plot(Omega_end) % plot final parameter set

% Set-Membership estimated system
dTheta_hat_end = dTheta_hat(end,:)';
AB_vec = vec(AB0);
AB_vec(idxJ) = AB_vec(idxJ) + dt * dTheta_hat_end;
AB_ms = reshape(AB_vec, size(AB0)); clear AB_vec dTheta_hat_end

A_ms = AB_ms(:,1:nx); B_ms = AB_ms(:,nx+1:nx+nu);
AB_ms = [A_ms B_ms]

% Output true system (discrete) for comparison
AB_true = [A_true B_true]

% Plot delta and absolute parameter values
for iP = 1:np
    figure(iP+10), %clf;
    subplot(1,2,1) % Delta parameter estimation
    hold on;
    patch([stepIdx'; flipud(stepIdx')], [dTheta_bounds(:,2*iP-1); flipud(dTheta_bounds(:,2*iP))], 'k', 'FaceAlpha', 0.1); % Prediction intervals
    plot(0:nSteps, dTheta_hat(:,iP), 'b');
    ylim([mean(dTheta_bounds(:,2*iP)), mean(dTheta_bounds(:,2*iP-1))])
    title(['Param ' num2str(iP) ' dTheta hat'])
    
    subplot(1,2,2) % Parameter estimation
    hold on
    % True value (linearized dynamics) dashed line
    plot(stepIdx, repmat(AB_true(idxJ(iP)), nSteps+1), 'k--');
    % Estimation history
    patch([stepIdx'; flipud(stepIdx')], AB0(idxJ(iP)) + dt * [dTheta_bounds(:,2*iP-1); flipud(dTheta_bounds(:,2*iP))], 'k', 'FaceAlpha', 0.1); % Prediction intervals
    plot(stepIdx, AB0(idxJ(iP)) + dt * dTheta_hat(:,iP), 'b');
    ylim(AB0(idxJ(iP)) + dt * [mean(dTheta_bounds(:,2*iP)), mean(dTheta_bounds(:,2*iP-1))])
    title(['Param ' num2str(iP) ' value'])
end

disp('Indices changed in A and B:')
disp([A_ms ~= A0, B_ms ~= B0])

sim.addModels(AB0, AB_ms);                          % add initial guess model and estimated model to simulator class
sim.simulate_estimSyst(t0, t_end, mdl.x_trim, U);   % simulate initial guess and estimated model
sim.plotEstimStateTrajectory();                     % plot state trajectories of true, initial and estimated system





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

