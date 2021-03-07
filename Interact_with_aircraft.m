function [x_ns] = Interact_with_aircraft(x_s,u_s,time_res)
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
    % Discretized system
    A = eye(4) + time_res*Ac;
    B = time_res*Bc;
    x_ns = A*x_s+B*[u_s];
end
