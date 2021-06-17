function [U,X,u] = MPC_tight(A,B,x_ini,ref,s,H,F,X_set, K)
    
    sx=4; %num. states
    su=1; %num. inputs
    x_i=sdpvar(sx, H+1); %state
    u_i=sdpvar(su, H); %input
    x_0=sdpvar(sx,1); %initial state
    
    %MPC
    objective_MPC=0;
    constraints_MPC=[x_i(:,1)==x_0];
    for i=1:H
        % Stage cost
        for it_s = 1:size(s,1)
            objective_MPC=objective_MPC   +   (x_i(s(it_s),i)-ref(it_s,i))^2;
        end
        % State Propagation Constraints
        constraints_MPC=[constraints_MPC, x_i(:,i+1) == A*x_i(:,i)+B*u_i(:,i)];
        
        % State & Input Constraints
%         constraints_MPC=[constraints_MPC, x_i(1,i+1) >= 0]; %Va (velocity norm) major to zero always
        % constraints_MPC=[constraints_MPC, deg2rad(-30) <= u_i(1,i) <= deg2rad(30)]; % Elevator deflection is limited by +-30 deg
        U_set = Polyhedron([1;-1], [deg2rad(50); deg2rad(50)]);     % Elevator deflection is limited by +-50 deg
        if(K==0)
            U_t{i} = U_set;
        else
            U_t{i} = U_set - K*F{i};
            U_t{i}.minHRep;
        end

        constraints_MPC=[constraints_MPC, U_t{i}.A * u_i(:,i) <= U_t{i}.b]; 
        X_t{i+1} = X_set - F{i+1};
        X_t{i+1}=X_t{i+1}.minHRep;
        constraints_MPC=[constraints_MPC, X_t{i+1}.A * x_i(:,i+1) <= X_t{i+1}.b];
    end
    ops = sdpsettings('verbose',1,'solver','sedumi');
    obj.optimizer=optimizer(constraints_MPC, objective_MPC, ops, {x_0}, {u_i,x_i});
    [sol, flag] = obj.optimizer(x_ini);
    U = sol{1};
    u=U(:,1);
    X = sol{2};
end
