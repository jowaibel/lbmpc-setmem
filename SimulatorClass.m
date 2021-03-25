classdef SimulatorClass < handle
    properties
        dt
        nx, nu
        dyn_func
        mdl
        A, B, x_trim, u_trim
        time, state_traj, lin_state_traj, lin_dyn_traj
    end
    
    methods
        function obj = SimulatorClass(nonlin_dyn_func, lin_mdl, dt)
            obj.mdl = lin_mdl;
            obj.A = obj.mdl.sys.A;
            obj.B = obj.mdl.sys.B;
            obj.x_trim = obj.mdl.x_trim; obj.nx = length(obj.x_trim);
            obj.u_trim = obj.mdl.u_trim; obj.nu = length(obj.u_trim);
            
            obj.dyn_func = nonlin_dyn_func;
            
            obj.dt = dt;
        end
        
        function [x] = simulate_one_step(obj, x0, u)
            % Sim nonlinear system one step
            [~, X_] = ode45( @(t_, x_) obj.dyn_func(x_, u), [0 obj.dt/2 obj.dt], x0);
            x = X_(end,:)';
        end
        function [xl, dxl] = simulate_one_step_lin(obj, x0, u)
            % Sim linear system one step
            [~, Xl_] = ode45( @(t_, x_) obj.A*(x_-obj.x_trim)+ obj.B*(u-obj.u_trim), [0 obj.dt/2 obj.dt], x0);
            xl = Xl_(end,:)';
            dxl = obj.A*(x0-obj.x_trim) + obj.B*(u-obj.u_trim);
        end
        
        function obj = simulate(obj, t0, t_end, x0, U)
            
            nSteps = ceil(t_end/obj.dt);
            
            T = (t0:obj.dt:t_end);
            X = [x0, zeros(4, nSteps)]; % nonlinear state trajectory
            Xl = X;                % linear state trajectory
            dXl = Xl;              % continuous linear dynamics trajectory
            
            for iStep = 1:nSteps
                X(:,iStep+1)  = obj.simulate_one_step(X(:,iStep), U(:,iStep));
                [Xl(:,iStep+1), dXl(:,iStep)] = obj.simulate_one_step_lin(Xl(:,iStep), U(:,iStep));
            end
            
            obj.time = T;
            obj.state_traj = X;
            obj.lin_state_traj = Xl;
            obj.lin_dyn_traj = dXl;
        end
        
        function plotStateTrajectrory(obj, T, X, Xl)
            
            if nargin < 2
                T = obj.time;
                X = obj.state_traj;
                Xl = obj.lin_state_traj;
            end
            
            figure(1);
            clf(1);
            for ix=1:obj.nx
                subplot(obj.nx, 1, ix); plot(T, X(ix,:), T, Xl(ix,:)); title(obj.mdl.sys.StateName{ix}), ylabel(obj.mdl.sys.StateUnit{ix});
            end
            legend('nonlin', 'lin');
        end
    end
end

