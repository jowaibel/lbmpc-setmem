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

% Check discretization uncertainty
A_true = sim.Ade; % Discretization from exact discretization
B_true = sim.Bde;
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

%% Initial guess for dynamics: XFLR model
JA = logical([
    1 1 1 1;
    1 1 1 1;
    1 1 1 0;
    0 0 1 0]);
JB = logical([
    0;
    0;
    1;
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

H_theta=[eye(np); -eye(np)];

% Bounds on delta parameters
%theta_rel_uncert = 1.0;
% Ac0Bc0 = [Ac0, Bc0];
% Ac0Bc0_pos = Ac0Bc0; Ac0Bc0_pos(Ac0Bc0 < 0) = 0;                % ### Asymmetric uncertainty
% Ac0Bc0_neg = Ac0Bc0; Ac0Bc0_neg(Ac0Bc0 > 0) = 0;                % ### Asymmetric uncertainty
% h_theta = theta_uncert * [Ac0Bc0_pos(idxJ); -Ac0Bc0_neg(idxJ)]; % ### Asymmetric uncertainty
%AB0_pos = AB0; AB0_pos(AB0 < 0) = 0;                % ### Asymmetric uncertainty (Discrete system)
%AB0_neg = AB0; AB0_neg(AB0 > 0) = 0;                % ### Asymmetric uncertainty
%h_theta = ( [-AB0_neg(idxJ); AB0_pos(idxJ)] + theta_rel_uncert * [AB0_pos(idxJ); -AB0_neg(idxJ)] ) ./ dt; % ### Bounds to avoid sign flip + uncertainty range in same sign
% h_theta = repmat(theta_uncert * abs(Ac0Bc0(idxJ)), 2, 1); % ### Symmetric uncertainty


% first option:
first_option = false;
theta_uncert = 1.0/dt;
AB0_parms=AB0(idxJ);
AB_true = [A_true B_true];
h_theta_upper=theta_uncert * ones(size(idxJ));
h_theta_lower=theta_uncert * ones(size(idxJ));
h_theta_upper(AB0_parms<0)=-AB0_parms(AB0_parms<0)/dt;
h_theta_lower(AB0_parms<0)=-AB0_parms(AB0_parms<0)/dt;

h_theta_lower(AB0_parms>0)=AB0_parms(AB0_parms>0)/dt;
h_theta_upper(AB0_parms>0)=AB0_parms(AB0_parms>0)/dt;
h_theta=[h_theta_upper;h_theta_lower];
header = {'lower bound','initial value','upper bound','true value','true value inside?'};
h_inf=[-h_theta_lower*dt+AB0(idxJ),AB0(idxJ),h_theta_upper*dt+AB0(idxJ),...,
    AB_true(idxJ),-h_theta_lower*dt+AB0(idxJ)<AB_true(idxJ) & AB_true(idxJ)<h_theta_upper*dt+AB0(idxJ)]; % to see params bound and see if AB_true is inside
output = [header; num2cell(h_inf)]; % run to see bounds

% second option:

theta_uncert = 1.0/dt;
AB0_parms=AB0(idxJ);
h_theta_upper=theta_uncert * ones(size(idxJ));
h_theta_lower=theta_uncert * ones(size(idxJ));
h_theta=[h_theta_upper;h_theta_lower];
header = {'lower bound','initial value','upper bound','true value','true value inside?'};
h_inf=[-h_theta_lower*dt+AB0(idxJ),AB0(idxJ),h_theta_upper*dt+AB0(idxJ),...,
    AB_true(idxJ),-h_theta_lower*dt+AB0(idxJ)<AB_true(idxJ) & AB_true(idxJ)<h_theta_upper*dt+AB0(idxJ)]; % to see params bound and see if AB_true is inside
output = [header; num2cell(h_inf)];  % run to see bounds

% Get Set Membership estimation
nSteps = nData;
nSteps = 100;
% Select nSteps samples from dataset, either randomly or equally distributed
%dataIdx = randperm(nData);          % Random selection from dataset
%dataIdx = 1:floor(nSteps/nSteps):nSteps; dataIdx = dataIdx(1:nSteps);  % Equally distributed selection over time (= downsampling)
dataIdx = 1:1:nSteps; % first nSteps points from dataset

nSteps = length(dataIdx);
stepIdx = 0:nSteps;

% Select data for estimation
DU = df_u(dataIdx, :)' - u_trim;        % control dataset
DX = df_lin_s(dataIdx, :)' - x_trim;    % linearized model, state dataset
DXP = df_lin_ns(dataIdx, :)' - x_trim;  % linearized model, next state dataset

AB_ms = AB0;
dTheta_hat = zeros(nSteps+1, np); % estimated values
dTheta_bounds = zeros(nSteps+1, 2*np);
setD = cell(1, nSteps+1);

% Define additive polytopic uncertainty description
DW = DXP - AB_ms * [DX; DU];
w_max = max(abs(DW), [], 2);
%w_max = 0.1 * ones(nx, 1); %0.041;
Hw = [eye(nx); -eye(nx)];

f1 = figure; iterations = 0; w_maxes = [];
clear dTheta_final_bounds_last

%%
recursive_estimation = false;
estimate_based_W = false;

term_crit = 10; % The estimation tries to tighten the dTheta uncertainty bounds until the certainty range in all parameters decreases less than term_crit.
        
if exist('dTheta_final_bounds_last', 'var'), fprintf('Warmstarting'); end

while (true)
    if first_option
        theta_uncert = 1.0/dt;
        AB0_parms=AB_ms(idxJ);
        h_theta_upper=theta_uncert * ones(size(idxJ));
        h_theta_lower=theta_uncert * ones(size(idxJ));
        h_theta_upper(AB0_parms<0)=-AB0_parms(AB0_parms<0)/dt;
        h_theta_lower(AB0_parms<0)=-AB0_parms(AB0_parms<0)/dt;

        h_theta_lower(AB0_parms>0)=AB0_parms(AB0_parms>0)/dt;
        h_theta_upper(AB0_parms>0)=AB0_parms(AB0_parms>0)/dt;
        h_theta=[h_theta_upper;h_theta_lower];
    end
    % define initial parameter set
    Omega{1} = Polyhedron(H_theta, h_theta);
    Hh_theta = H_theta .* h_theta;
    dTheta_bounds(1,:) = Hh_theta(logical(repmat(eye(np),2,1))); clear Hh_theta
    
    % Define disturbance bounds
    if estimate_based_W
        DW = DXP - AB_ms * [DX; DU];
        w_max = max(abs(DW), [], 2);
    end
    hw = repmat(w_max, 2, 1);
    W = Polyhedron(Hw, hw);
    
    % Instantiate set membership estimator
    sm = SetMembership(Omega{1}, W, ABi, AB0);
    
    % Keep track of w_max evolution
    w_maxes = [w_maxes w_max]; 
    if size(w_maxes,2)-1 > iterations(end), iterations = [iterations iterations(end)+1]; end
    set(0, 'currentfigure', f1);
    subplot(length(w_max), 1, 1); plot(iterations, w_maxes(1,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(iterations), grid on, title('w\_{max}')
    subplot(length(w_max), 1, 2); plot(iterations, w_maxes(2,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(iterations), grid on
    subplot(length(w_max), 1, 3); plot(iterations, w_maxes(3,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(iterations), grid on
    subplot(length(w_max), 1, 4); plot(iterations, w_maxes(4,:), '.-k', 'MarkerSize', 12), xlabel('it'), xticks(iterations), grid on, drawnow

    fprintf(['\nStarting SM estimation with w_max = ' num2str(w_max')]); tic;    
    if recursive_estimation
        % Estimate recursively, sample by sample
        
        % Use nSteps samples for set membership identification
        iPlot = 1;
        for iStep = stepIdx(2:end)
            id = dataIdx(iStep);
            [Omega{iStep+1}, setD{iStep+1}] = sm.update(DXP(:,id), DX(:,id), DU(:,id));    % update set membership
            
            if Omega{iStep+1}.isEmptySet
                if estimate_based_W
                    AB_ms = sm.get_AB(); % Save current AB estimate
                    modify_W_str = 'iterate';
                else
                    % Restart estimation with larger W
                    w_max = 1.1 * w_max;
                    modify_W_str = 'enlarge';
                end
                fprintf([' -- Step ' num2str(iStep) '/' num2str(nSteps) ': Empty set, ' modify_W_str ' disturbance set W.']);
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
        Omega{nSteps+1} = sm.update(DXP, DX, DU);
        if Omega{nSteps+1}.isEmptySet
            if estimate_based_W
                AB_ms = sm.get_AB(); % Save current AB estimate
                modify_W_str = 'iterate';
            else
                % Restart estimation with shrinked W
                w_max = 1.1 * w_max;
                modify_W_str = 'enlarge';
            end
            fprintf(['\nAll samples used. Empty set, ' modify_W_str ' disturbance set W.']);
        else
           
            dTheta_hat(nSteps+1,:) = sm.theta_hat';  % estimate parameter (center of the estimated set)
            dTheta_hat(:,:) = interp1([1 2],[dTheta_hat(1,:);dTheta_hat(end,:)], linspace(1,2,nSteps+1));

            dTheta_bounds(nSteps+1,:) = sm.theta_bounds';
            dTheta_bounds(:,:) = interp1([1 2],[dTheta_bounds(1,:);dTheta_bounds(end,:)], linspace(1,2,nSteps+1));
            iStep = nSteps;
        end
    end
    
    if iStep == nSteps
        AB_ms = sm.get_AB(); % Save current AB estimate
        % Prepare restart of estimation with updated W (warmstart)
        if estimate_based_W
            modify_W_str = 'Iterate';
        else
            % Restart estimation with shrinked W
            w_max = 0.9 * w_max;
            modify_W_str = 'Tighten';
        end
        
        if exist('dTheta_final_bounds_last', 'var')
            % Range between dTheta bounds
            dTheta_final_bounds_size = abs(diff(reshape(dTheta_bounds(end,:), 2, np)));
            dTheta_final_bounds_last_size = abs(diff(reshape(dTheta_final_bounds_last, 2, np)));
            
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
        fprintf(['\nAll samples used (' num2str(round(toc,1)) 's). ' modify_W_str ' disturbance set W.\n ----']);
        dTheta_final_bounds_last = dTheta_bounds(end,:);
    end
end
%Omega_end = Omega{nSteps+1}; %figure, plot(Omega_end) % plot final parameter set

% Set-Membership estimated system
AB_ms
A_ms = AB_ms(:,1:nx); B_ms = AB_ms(:,nx+1:nx+nu);

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





%return

% These membership estimations should'nt converge, as the initial guess is the horizontal line, its a different plot that the one used before.



% figure
% sgtitle("Nextstate v/s state, Linearized model (red), data (blue), M.S model (green)")
% states=["Va","alpha","q","gamma"];
% for var = 1:4
%     subplot(2,2,var)
%     plot(df_s(:,var),df_ns(:,var), 'bp', 'MarkerSize', 0.5) %state against next state
%     hold on
%     plot(df_s(:,var),df_lin_ns(:,var), 'rp', 'MarkerSize', 0.5)
%     hold on
%     plot(df_s(:,var),df_ms_ns(:,var), 'gp', 'MarkerSize', 1)
%     xlim([min(df_s(:,var)),max(df_s(:,var))]);
%     ylim([min(df_ns(:,var)),max(df_ns(:,var))]);
%     xlabel(states(var));
% end
% 
% figure
% sgtitle(["Nextstate(data) - Nextstate(linearized model) [blue]"; ...
%     "and Nextstate(data) - Nextstate(membership model) [green]"])
% for var = 1:4
%     subplot(2,2,var)
%     plot(df_s(:,var),df_ns(:,var)-df_lin_ns(:,var), 'bp', 'MarkerSize', 1)
%     hold on;
%     plot(df_s(:,var),df_ns(:,var)-df_ms_ns(:,var), 'gp', 'MarkerSize', 1)
%     xlim([min(df_s(:,var)),max(df_s(:,var))]);
%     xlabel(states(var));
% end

%% Maneuver with linearized model and membership model the non linear simulator
see_progress=true; %show progress bar
AB_ms=[A_ms B_ms];

% Set up reference
ref_q = [repmat(deg2rad(0),1,40)];
ref_theta = [repmat(deg2rad(5),1,100)];
ref=[ref_theta];
s=[4]; %state(s) to control in concordance with ref order (in this case applying only theta control)
H=5; % MPC horizon
opt_steps=size(ref,2)-H;
disp("Controlling state variable ["+char(mdls.xflr_uw.sys.StateName(s))+"] for "+dt*size(ref,2)+" seconds, horizon H="+H+" steps ("+dt*H+" seconds)")

% Calculate noise bounds:
DW_ms = DXP - AB_ms * [DX; DU];
DW_true = DXP - AB_true * [DX; DU];
w_max_ms = max(abs(DW_ms), [], 2);
w_max_true = max(abs(DW_true), [], 2);
Hw = [eye(nx);-eye(nx)];
W_ms = Polyhedron(Hw,[w_max_ms; w_max_ms]);
W_true = Polyhedron(Hw,[w_max_true; w_max_true]);

% Ancillary/Tube Controller for Constraint Tightening
use_tube_controller = false;
if use_tube_controller
    [K_ms, P_ms] = dlqr(A_ms,B_ms,eye(nx),10*eye(nu));     % Q = 10*eye(nx), R = 10*eye(nu)
    K_ms = -K_ms;
    [K_true, P_true] = dlqr(A_true,B_true,eye(nx),10*eye(nu));     % Q = 10*eye(nx), R = 10*eye(nu)
    K_true = -K_true;
else
    K_ms = 0;
    K_true = 0;
end

% Calculate disturbance reachable sets
upper_ms=[0;0;0;0];
lower_ms=[0;0;0;0];
upper_true=[0;0;0;0];
lower_true=[0;0;0;0];
F_ms{1} = Polyhedron(zeros(0,nx),[]);   %start with empty poylhedron
F_true{1} = Polyhedron(zeros(0,nx),[]);   %start with empty poylhedron
if see_progress progressbar= waitbar(0, 'Starting'); end % Init single bar
for i=1:H
    if use_tube_controller
        F_ms{i+1} = plus((A_ms+B_ms*K_ms)*F_ms{i}, W_ms, 'vrep');
        F_true{i+1} = plus((A_true+B_true*K_true)*F_true{i}, W_true, 'vrep');
    else
        F_ms{i+1} = plus(A_ms*F_ms{i},W_ms,'vrep');
        F_true{i+1} = plus(A_true*F_true{i},W_true,'vrep');
    end
    F_ms{i+1}.minHRep;
    F_true{i+1}.minHRep;
    if i>1
        proj_ms=[projection(F_ms{i},[1]).H ; projection(F_ms{i},[2]).H ; projection(F_ms{i},[3]).H ; projection(F_ms{i},[4]).H];
        upper_ms=[upper_ms, abs(proj_ms(1:2:end,2)./proj_ms(1:2:end,1))];
        lower_ms=[lower_ms, -abs(proj_ms(2:2:end,2)./proj_ms(2:2:end,1))];
        
%         proj_true=[projection(F_true{i},[1]).H ; projection(F_true{i},[2]).H ; projection(F_true{i},[3]).H ; projection(F_true{i},[4]).H];
%         upper_true=[upper_true, abs(proj_true(1:2:end,2)./proj_true(1:2:end,1))];
%         lower_true=[lower_true, -abs(proj_true(2:2:end,2)./proj_true(2:2:end,1))];
    end
    if see_progress
        waitbar(i/H, progressbar, sprintf('Dist. reachable set Progress: %d %%', floor(i/H*100))); % Update progress bar
    end
end
if see_progress close(progressbar); end

figure
for k =1:nx
    subplot(4,2,k*2-1);fill([1:H, fliplr(1:H)],[lower_ms(k,:), fliplr(upper_ms(k,:))],'g');
    if k==1 title("Disturb. reach. sets M.S"); end
    ylabel(char(mdls.xflr_uw.sys.StateName(k)));
end
% for k =1:nx
%     subplot(4,2,k*2);fill([1:H, fliplr(1:H)],[lower_true(k,:), fliplr(upper_true(k,:))],'g');
%     if k==1 title("Disturb. reach. sets true"); end
%     ylabel(char(mdls.xflr_uw.sys.StateName(k)));
% end
%% 
% specify states domain
H_x=[eye(nx);-eye(nx)];
% h_x=[13.5 ; 0 ; .4 ; .2; 
%      -12 ; 3 ; .4 ; .2];
h_x = [25; 5; 10; deg2rad(70);       
       -5; 5; 10; deg2rad(70)]; 
X_set=Polyhedron(H_x,h_x);

%with membership model
x0=x_trim;
X_ms=x0;
if see_progress progressbar= waitbar(0, 'Starting'); end
for t = 1:opt_steps
    [U,X,u]=MPC_tight(A_ms,B_ms,x0,ref(:,t:end),s,H,F_ms,X_set,K_ms);
%     x0=A_true*x0+B_true*u;
    x0 = sim.simulate_one_step(x0, u);
    X_ms=cat(2,X_ms,x0);
    if isnan(u)
        fprintf('\nMPC with est. model infeasible');
        break
    end
    if see_progress
        waitbar(t/opt_steps, progressbar, sprintf('M.S model MPC Progress: %d %%', floor(t/opt_steps*100))); % Update progress bar
    end
end
if see_progress close(progressbar); end

%with linearized model
x0=x_trim;
X_true=x0;
if see_progress progressbar= waitbar(0, 'Starting'); end % Init single bar
for t = 1:opt_steps
    [U,X,u]=MPC_tight(A_true,B_true,x0,ref(:,t:end),s,H,F_true,X_set,K_true);
%     x0=A_true*x0+B_true*u;
    x0 = sim.simulate_one_step(x0, u);
    X_true=cat(2,X_true,x0);
    if isnan(u)
        fprintf('\nMPC with true model infeasible');
        break
    end
    if see_progress
        waitbar(t/opt_steps, progressbar, sprintf('True model MPC Progress: %d %%', floor(t/opt_steps*100))); % Update progress bar
    end
end
if see_progress close(progressbar); end

figure
for k=1:nx
    subplot(4,2,2*k-1);plot(X_ms(k,:));hold on;ylabel(char(mdls.xflr_uw.sys.StateName(k)));
    if k==1 title("MPC M.S constr. satisf."); end
    plot(1:opt_steps,h_x(k)*ones(opt_steps,1),"--");hold on;plot(1:opt_steps,-h_x(k+4)*ones(opt_steps,1),"--")
end
for k=1:nx
    subplot(4,2,2*k);plot(X_true(k,:));hold on;ylabel(char(mdls.xflr_uw.sys.StateName(k)));
    if k==1 title("MPC true model constr. satisf."); end
    plot(1:opt_steps,h_x(k)*ones(opt_steps,1),"--");hold on;plot(1:opt_steps,-h_x(k+4)*ones(opt_steps,1),"--")
end

figure
for it_s = 1:size(s,1)
    subplot(2+size(s,1),1,it_s)
    plot(X_true(s(it_s),1:end),"r")
    hold on
    plot(X_ms(s(it_s),1:end),"g")
    hold on
    plot(ref(it_s,:),"k")
    legend(["True","M.S","ref."])
    ylabel(char(mdls.xflr_uw.sys.StateName(s(it_s))))
end

obj_true=sum((X_true(s,:)-ref(:,1:size(X_true,2))).^2,1);
obj_ms=sum((X_ms(s,:)-ref(:,1:size(X_ms,2))).^2,1);
subplot(2+size(s,1),1,size(s,1)+1)
plot(obj_true,"r")
hold on
plot(obj_ms,"g")
legend(["True","M.S"])
title("Objective function evaluation")

subplot(2+size(s,1),1,size(s,1)+2)
cs_true=cumsum(obj_true);
cs_ms=cumsum(obj_ms);
plot(cs_true,"r")
hold on
plot(cs_ms,"g")
legend(["True: "+cs_true(end),"M.S: "+cs_ms(end)])
title("Cumulative objective function evaluation")
