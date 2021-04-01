classdef SimulatorClass < handle
    properties
        dt
        nx, nu
        dyn_func
        mdl
        Ac, Bc, x_trim, u_trim, Ad, Bd, Ade, Bde,
        time, state_traj, lin_state_traj, lin_dyn_traj, lind_state_traj, linde_state_traj
    end
    
    methods
        function obj = SimulatorClass(nonlin_dyn_func, lin_mdl, dt)
            obj.dt = dt;
            obj.x_trim = lin_mdl.x_trim; obj.nx = length(obj.x_trim);
            obj.u_trim = lin_mdl.u_trim; obj.nu = length(obj.u_trim);
            
            % Nonlinear model
            obj.dyn_func = nonlin_dyn_func;
            
            % Linearized model
            obj.mdl = lin_mdl;
            obj.Ac = obj.mdl.sys.A;
            obj.Bc = obj.mdl.sys.B;
            
            % (Euler) >D<iscretization -> d
            obj.Ad = eye(obj.nx) + obj.Ac * dt; % Discretization (Euler)
            obj.Bd = obj.Bc * dt; % Discretization (Euler)
            
            % >E<xact >D<iscretization -> de
            M = [obj.Ac obj.Bc; zeros(obj.nu, obj.nx+obj.nu)];
            F = expm(M*dt);
            obj.Ade = F(1:obj.nx, 1:obj.nx);     % Exact discretization
            obj.Bde = F(1:obj.nx, obj.nx+1:end); % Exact discretization
        end
        
        function [x] = simulate_one_step(obj, x0, u)
            % Sim nonlinear system one step
            [~, X_] = ode45( @(t_, x_) obj.dyn_func(x_, u), [0 obj.dt/2 obj.dt], x0);
            x = X_(end,:)';
        end
        function [xl, dxl] = simulate_one_step_lin(obj, x0, u)
            % Sim linear system one step
            [~, Xl_] = ode45( @(t_, x_) obj.Ac*(x_-obj.x_trim)+ obj.Bc*(u-obj.u_trim), [0 obj.dt/2 obj.dt], x0);
            xl = Xl_(end,:)';
            dxl = obj.Ac*(x0-obj.x_trim) + obj.Bc*(u-obj.u_trim);
        end
        function [xp] = simulate_one_step_lind(obj, x0, u)
            % Sim linear discrete system one step
            xp = obj.x_trim + obj.Ad * (x0-obj.x_trim) + obj.Bd * (u-obj.u_trim);
        end
        function [xp] = simulate_one_step_linde(obj, x0, u)
            % Sim linear discrete (exact) system one step
            xp = obj.x_trim + obj.Ade * (x0-obj.x_trim) + obj.Bde * (u-obj.u_trim);
        end
        
        function obj = simulate(obj, t0, t_end, x0, U)
            
            nSteps = ceil(t_end/obj.dt);
            
            T = (t0:obj.dt:t_end);
            X = [x0, zeros(4, nSteps)]; % nonlinear state trajectory
            Xl = X;                % linear state trajectory
            dXl = Xl;              % continuous linear dynamics trajectory
            Xld = [x0, zeros(4, nSteps)]; % linear discrete state trajectory
            Xlde = [x0, zeros(4, nSteps)]; % linear discrete state trajectory
            
            for iStep = 1:nSteps
                X(:,iStep+1)  = obj.simulate_one_step(X(:,iStep), U(:,iStep));
                [Xl(:,iStep+1), dXl(:,iStep)] = obj.simulate_one_step_lin(Xl(:,iStep), U(:,iStep));
                Xld(:,iStep+1) = obj.simulate_one_step_lind(Xld(:,iStep), U(:,iStep));
                Xlde(:,iStep+1) = obj.simulate_one_step_linde(Xlde(:,iStep), U(:,iStep));
            end
            
            obj.time = T;
            obj.state_traj = X;
            obj.lin_state_traj = Xl;
            obj.lin_dyn_traj = dXl;
            obj.lind_state_traj = Xld;
            obj.linde_state_traj = Xlde;
        end
        
        function plotStateTrajectrory(obj, T, X, Xl)
            
            if nargin < 2
                T = obj.time;
                X = obj.state_traj;
                Xl = obj.lin_state_traj;
                Xld = obj.lind_state_traj;
                Xlde = obj.linde_state_traj;
            end
            
            figure(1);
            clf(1);
            for ix=1:obj.nx
                subplot(obj.nx, 1, ix);
                plot(T, X(ix,:), T, Xl(ix,:), T, Xld(ix,:), '.', T, Xlde(ix,:), '.');
                title(obj.mdl.sys.StateName{ix});
                ylabel(obj.mdl.sys.StateUnit{ix});
            end
            subplot(obj.nx, 1, 1);
            legend('nonlin', 'lin', 'lind', 'linde');
        end
    end
end

