function [U,X,u] = MPC_nom_lin(x_ini,time_res,alpha_ref,gamma_ref,H)
    Ac = [ 
      -0.147494    -3.90244           0    -9.76268;
      -0.123791     -9.8892           1   0.0733121;
    0.000229919    -94.8597    -19.1447           0;
       0.123791      9.8892           0  -0.0733121];
    Bc = [ 
               0;
               0;
        -103.768;
               0];
    xtrim=[12.5591, -0.0721008, -1.31188e-08, -0.0940332]';
    utrim=-0.0110968;
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
        objective_MPC=objective_MPC   +   (x_i(2,i+1)-alpha_ref(1,i))^2   +    (x_i(4,i+1)-gamma_ref(1,i))^2;
        % State Propagation Constraints
        constraints_MPC=[constraints_MPC, x_i(:,i+1) == time_res*(Ac*(x_i(:,i)-xtrim)+Bc*(u_i(:,i)-utrim))+x_i(:,i)];
        
        % State & Input Constraints
        constraints_MPC=[constraints_MPC, x_i(1,i+1) >= 0]; %Va (velocity norm) major to zero always
        constraints_MPC=[constraints_MPC, deg2rad(-30) <= u_i(1,i) <= deg2rad(30)]; % Elevator deflection is limited by +-30 deg
    end
    ops = sdpsettings('verbose',1,'solver','sedumi');
    obj.optimizer=optimizer(constraints_MPC, objective_MPC, ops, {x_0}, {u_i,x_i});
    [sol, flag] = obj.optimizer(x_ini);
    U = sol{1};
    u=U(:,1);
    X = sol{2};
end
