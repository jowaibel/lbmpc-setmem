classdef SimulatorClass < handle
    properties
        dt
        nx, nu
        U
        dyn_func
        mdl
        x_trim, u_trim
        time
        Ac, Bc              % continuous system matrices
        Ad, Bd              % discrete system matrices Ade, Bde
        Ade, Bde
        nonlin_state_traj   % state trajectory from nonlinear continuous simulation
        linc_state_traj     % state trajectory from linear continuous simulation
        linc_dyn_traj       % dynamics trajectory from linear continuous simulation
        lind_state_traj     % state trajectory from linear discrete simulation
        linde_state_traj    % state trajectory from linear discrete (exact) simulation
        process_noise_abs
    end
    
    methods
        function obj = SimulatorClass(nonlin_dyn_func, lin_mdl, dt, process_noise_abs)
            
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
            
            if nargin < 4
                process_noise_abs = zeros(obj.nx, 1);
            end
            obj.process_noise_abs = process_noise_abs;
        end

        % Steppers for nonlinear / linear cont. / linear discr, dynamics
        function [x] = simulate_nonlinear_step(obj, dyn_func, x0, u, dt)
            % Sim nonlinear system for dt
            [~, X_] = ode45( @(t_, x_) dyn_func(x_, u), [0 dt/2 dt], x0);
            x = X_(end,:)';
        end
        
        function [xl, dxl] = simulate_lin_cont_step(obj, Ac, Bc, x_trim, u_trim, x0, u, dt)
            % Sim linear system for dt
            dyn_func = @(t_, x_) Ac * (x_ - x_trim) + Bc * (u - u_trim);

            [~, Xl_] = ode45(dyn_func, [0 dt/2 dt], x0);
            xl = Xl_(end,:)';
            dxl = dyn_func(dt, x0);

        end
        
        function [xp] = simulate_lin_discr_step(obj, Ad, Bd, x_trim, u_trim, x0, u)
            % Sim linear discrete system
            xp = x_trim + Ad * (x0 - x_trim) + Bd * (u - u_trim);
        end
        
        % Integration along input trajectory
        function [X] = simulate_nonlinear(obj, dyn_func, dt, x0, U, D)
            
            nSteps = size(U, 2);
            X = [x0, zeros(size(x0, 1), nSteps)];
            
            for iStep = 1:nSteps
                X(:, iStep+1)  = obj.simulate_nonlinear_step(dyn_func, X(:,iStep), U(:,iStep), dt) + D(:,iStep);
            end
        end
        
        function [X, DX] = simulate_linear_cont(obj, Ac, Bc, x_trim, u_trim, dt, x0, U, D)
            
            nSteps = size(U, 2);
            X = [x0, zeros(size(x0, 1), nSteps)];
            DX = [x0, zeros(size(x0, 1), nSteps)];
            
            for iStep = 1:nSteps
                [X(:, iStep+1), DX(:,iStep)] = obj.simulate_lin_cont_step(Ac, Bc, x_trim, u_trim, X(:,iStep), U(:,iStep), dt);
                X(:, iStep+1) = X(:, iStep+1) + D(:,iStep);
            end
        end
        
        function [X] = simulate_linear_discr(obj, Ad, Bd, x_trim, u_trim, x0, U, D)
            
            nSteps = size(U, 2);
            X = [x0, zeros(size(x0, 1), nSteps)];
            if nargin < 8
                D = zeros(size(X));
            end
            
            for iStep = 1:nSteps
                X(:, iStep+1)  = obj.simulate_lin_discr_step(Ad, Bd, x_trim, u_trim, X(:,iStep), U(:,iStep)) + D(:,iStep);
            end
        end
        
        % Simulation with noise along input trajectory and for different configurations
        function obj = simulate_all(obj, T, x0, U)
            
            nSteps = length(T) - 1;
            dt = T(2) - T(1);
            
            % Disturbance realisation along trajectory
            D = diag(obj.process_noise_abs) * ( 2 * (rand(obj.nx, length(T)) - 0.5) );
            obj.U=U;
            obj.time = T;
            
            obj.nonlin_state_traj = obj.simulate_nonlinear(obj.dyn_func, dt, x0, U(:, 1:end-1), D);
            
            [Xlinc, DXlinc] = obj.simulate_linear_cont(obj.Ac, obj.Bc, obj.x_trim, obj.u_trim, dt, x0, U(:, 1:end-1), D);
            obj.linc_state_traj = Xlinc;
            obj.linc_dyn_traj = DXlinc;
            
            obj.lind_state_traj = obj.simulate_linear_discr(obj.Ad, obj.Bd, obj.x_trim, obj.u_trim, x0, U(:, 1:end-1), D);
            obj.linde_state_traj = obj.simulate_linear_discr(obj.Ade, obj.Bde, obj.x_trim, obj.u_trim, x0, U(:, 1:end-1), D);
        end
        
        function plotStateTrajectory(obj, X_data)
            
            T = obj.time;
            X = obj.nonlin_state_traj;
            
            if nargin < 2
                % Plot nonlinear simulation and various linear/discretized
                Xlc = obj.linc_state_traj;
                Xld = obj.lind_state_traj;
                Xlde = obj.linde_state_traj;
                
                figure;
                for ix=1:obj.nx
                    subplot(obj.nx/2+1, 2, ix);
                    plot(T, X(ix,:), ...
                        T, Xlc(ix,:), ...
                        T, Xld(ix,:), '.', ...
                        T, Xlde(ix,:), '.');
                    xlim([0,T(end)]);
                    title(obj.mdl.sys.StateName{ix});
                    ylabel(obj.mdl.sys.StateUnit{ix});
                end
                subplot(obj.nx/2+1, 2, 1);
                legend('nonlin', 'lin', 'lind', 'linde');
                xlim([0,T(end)]);
            else
                % Plot data and nonlinear/linear simulation
                Xlc = obj.linc_state_traj;
                figure;
                
                for ix=1:obj.nx
                    subplot(obj.nx/2+1, 2, ix);
                    plot(T, X_data(ix,:), ...
                        T, X(ix,:), ...
                        T, Xlc(ix,:));
                    xlim([0,T(end)]);
                    title(obj.mdl.sys.StateName{ix});
                    ylabel(obj.mdl.sys.StateUnit{ix});
                end
                subplot(obj.nx/2+1, 2,1);
                xlim([0,T(end)]);
                legend('data', 'nonlin', 'lin');
            end
            subplot(obj.nx+1, 1, 5);
            plot(T,obj.U,"k")
            title("\delta_e");
            ylabel("rad");
            xlim([0,T(end)]);
            xlabel("Time (s)");
            subplot(obj.nx/2+1, 2, 3);xlabel("Time (s)");
            subplot(obj.nx/2+1, 2, 4);xlabel("Time (s)");
        end
        function plotStateTrajectoryComparison(obj, ...
                T1, X1, X1_legend, ...
                T2, X2, X2_legend, ...
                T3, X3, X3_legend, ...
                T4, X4, X4_legend)
            
            % Plot data and nonlinear/linear simulation
            figure;
            
            for ix = 1:obj.nx
                subplot(obj.nx, 1, ix);
                hold on
                
                plot(T1, X1(ix,:));
                
                if nargin > 4
                    plot(T2, X2(ix,:));
                end
                if nargin > 7
                    plot(T3, X3(ix,:));
                end
                if nargin > 10
                    plot(T4, X4(ix,:));
                end
                
                title(obj.mdl.sys.StateName{ix});
                ylabel(obj.mdl.sys.StateUnit{ix});
            end
            
            subplot(obj.nx, 1, 1);
            if nargin > 10
                legend(X1_legend, X2_legend, X3_legend, X4_legend);
            elseif nargin > 7
                legend(X1_legend, X2_legend, X3_legend);
            elseif nargin > 4
                legend(X1_legend, X2_legend);
            else
                legend(X1_legend);
            end
        end
        
        function [NRMSE] = calculateNRMSE(obj, X_data, X)
            RMSE = sqrt(sum((X_data-X).^2, 2));     % nx RMSE-values of the initial guess model
            NRMSE = mean(RMSE./abs(mean(X,2)));     % Mean of nx normalized RMSE-values
        end
    end
end