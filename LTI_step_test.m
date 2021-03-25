mdls = lon_LTI_models();

%% Induce step
sys_ss = mdls.uw.sys;

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