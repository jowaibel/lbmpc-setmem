addpath('nl_dynamics')

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
       
% Nonlinear simulation exmaple for 'gamma model'
t0 = 0;
t_end = 10;
dt = 0.002;
nSteps = ceil(t_end/dt);
xtrim = [12.5591, -0.0721008, -1.31188e-08, -0.0940332]'; % [Va alpha q gamma]
u_trim = [-0.0110968, 0]';
x0 = xtrim;

T = (t0:dt:t_end)';
X = [x0'; ...
    zeros(nSteps, 4)];
Xl = X;
dXl = Xl;

for iStep = 1:nSteps
    u = [-0.0110968 0];
    [~, X_] = ode45(@(t_, x_) dyn_func_gam(x_, u), [0 dt/2 dt], X(iStep,:));
    X(iStep+1,:) = X_(end,:);
    [~, Xl_] = ode45(@(t_, x_) A*(x_ - xtrim)+B*(u'- u_trim), [0 dt/2 dt], Xl(iStep,:)');
    Xl(iStep+1,:) = Xl_(end,:);
    dXl(iStep,:) = A*(Xl(iStep,:)'- xtrim)+B*u';
end


%%
figure
subplot(4, 1, 1); plot(T, X(:,1), T, Xl(:,1)); title('x_1')
subplot(4, 1, 2); plot(T, X(:,2), T, Xl(:,2)); title('x_2')
subplot(4, 1, 3); plot(T, X(:,3), T, Xl(:,3)); title('x_3')
subplot(4, 1, 4); plot(T, X(:,4), T, Xl(:,4)); title('x_4')
legend('nonlin', 'linzd');




