classdef SimulatorClass < handle
    properties
        dt
        nx, nu
        dyn_func
        mdl
        Ac, Bc, x_trim, u_trim, Ad, Bd, Ade, Bde, A0, B0, A_ms, B_ms
        time, state_traj, lin_state_traj, lin_dyn_traj, lind_state_traj, linde_state_traj
        lin_init_state_traj, lin_ms_state_traj
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
        
        function obj = addModels(obj, AB0, AB_ms)
           obj.A0 = AB0(:,1:obj.nx);                % initial guess (discrete)
           obj.B0 = AB0(:,obj.nx+1:obj.nx+obj.nu);
           obj.A_ms = AB_ms(:,1:obj.nx);              % estimated system (discrete)
           obj.B_ms = AB_ms(:,obj.nx+1:obj.nx+obj.nu);
            
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
        function [xp] = simulate_one_step_init(obj, x0, u)
            % Sim linear discrete initial guess system one step
            xp = obj.x_trim + obj.A0 * (x0-obj.x_trim) + obj.B0 * (u-obj.u_trim);
        end
        function [xp] = simulate_one_step_estim(obj, x0, u)
            % Sim linear discrete estimated system one step
            xp = obj.x_trim + obj.A_ms * (x0-obj.x_trim) + obj.B_ms * (u-obj.u_trim);
        end
        
        function obj = simulate(obj, t0, t_end, x0, U)
            
            nSteps = ceil((t_end-t0)/obj.dt);
            
            T = (t0:obj.dt:t_end);
            X = [x0, zeros(obj.nx, nSteps)]; % nonlinear state trajectory
            Xl = X;                % linear state trajectory
            dXl = Xl;              % continuous linear dynamics trajectory
            Xld = [x0, zeros(obj.nx, nSteps)]; % linear discrete state trajectory
            Xlde = [x0, zeros(obj.nx, nSteps)]; % linear discrete state trajectory
            
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
        
        function obj = simulate_estimSyst(obj, t0, t_end, x0, U)
            
            nSteps = ceil((t_end-t0)/obj.dt);
            
            Xld0 = [x0, zeros(obj.nx, nSteps)]; % linear discrete state trajectory (initial guess system)
            Xld_ms = [x0, zeros(obj.nx, nSteps)]; % linear discrete state trajectory (estimated system)
            
            for iStep = 1:nSteps
                Xld0(:,iStep+1) = obj.simulate_one_step_init(Xld0(:,iStep), U(:,iStep));
                Xld_ms(:,iStep+1) = obj.simulate_one_step_estim(Xld_ms(:,iStep), U(:,iStep));
            end
            
            obj.lin_init_state_traj = Xld0;     
            obj.lin_ms_state_traj = Xld_ms;
        end
        
        function plotStateTrajectory(obj, T, X, Xl)
            
            if nargin < 2
                T = obj.time;
                X = obj.state_traj;
                Xl = obj.lin_state_traj;
                Xld = obj.lind_state_traj;
                Xlde = obj.linde_state_traj;
            end
            
            figure;
            for ix=1:obj.nx
                subplot(obj.nx, 1, ix);
                plot(T, X(ix,:), T, Xl(ix,:), T, Xld(ix,:), '.', T, Xlde(ix,:), '.');
                title(obj.mdl.sys.StateName{ix});
                ylabel(obj.mdl.sys.StateUnit{ix});
            end
            subplot(obj.nx, 1, 1);
            legend('nonlin', 'lin', 'lind', 'linde');
        end
        
        function plotEstimStateTrajectory(obj, T, X0, X_ms)
            
            if nargin < 2
                T = obj.time;
                X0 = obj.lin_init_state_traj;
                X_ms = obj.lin_ms_state_traj;
                Xl = obj.lin_state_traj;
            end
            
            figure;
            for ix=1:obj.nx
                subplot(obj.nx, 1, ix);
                plot(T, Xl(ix,:), T, X0(ix,:),'--', T, X_ms(ix,:),'.-');
                title(obj.mdl.sys.StateName{ix});
                ylabel(obj.mdl.sys.StateUnit{ix});
            end
            subplot(obj.nx, 1, 1);
            legend('true','initial guess', 'estimated system');
        end
    end
end

