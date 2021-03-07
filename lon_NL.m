addpath('nl_dynamics')

% Nonlinear simulation exmaple for 'gamma model'
t0 = 0;
t_end = 10;
dt = 0.02;
nSteps = ceil(t_end/dt);
x0 = [12; 0; 0; 0]; % [Va alpha q gamma]

T = (t0:dt:t_end)';
X = [x0'; ...
    zeros(nSteps, 4)];

for iStep = 1:nSteps
    u = 0;
    func_wr = @(t_, x_) dyn_func_gam(x_, u);
    [~, X_] = ode45(func_wr, [0 dt/2 dt], X(iStep,:));
    X(iStep+1,:) = X_(end,:);
end
