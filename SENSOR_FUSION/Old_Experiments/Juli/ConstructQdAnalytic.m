%% Construct Qd based on analytical results from Sebastian Hilsenbeck 

function Qd= ConstructQdAnalytic(delta_t,C_q_w_i_hat, Cross_a_x, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)

    % Computediscrete -  time covariance matrix (Qd = Integrate[Fd (t)*Gc*Qc*Gc'*Fd' (t) dt])
    % var_na, var_nomega, etc.are noise varainces

    Qd = zeros (22);
    % left no scalar multiplying
    % all the dt^x scalars are multiplied by eye(3)
    % all the variances are diagonalized and moved to the correct? side of matrix multiplication 
    % C_q_w_i_hat'*Cross_a_x on the left, Cross_a_x '* C_q_w_i_hat' on the right
    % Don't we have to traspose Cross_a_x when it is at the right side? 

    dt =  delta_t;
    
    Qd (4 : 6,   1 : 3) =                             diag(var_na      *dt^2/2 ...
                                                    +      var_nba     *dt^4/8) ...
                            - C_q_w_i_hat'*Cross_a_x*((dt^4/8*eye(3) - dt^5/60*Cross_omega_x + dt^6/144*Cross_omega_x2) * diag(var_nomega) ...
                                                    + diag(var_nbomega) *(-dt^6/72*eye(3) + dt^7/1008*Cross_omega_x - dt^8/1920*Cross_omega_x2))*Cross_a_x*C_q_w_i_hat;
    %corrected
    Qd (7 : 9,   1 : 3) =     ((dt^3/6*eye(3) - dt^4/12*Cross_omega_x + dt^5/40*Cross_omega_x2)* diag(var_nomega) ...
                            + diag(var_nbomega)* (dt^5/30*eye(3) - dt^6/144*Cross_omega_x + dt^7/840*Cross_omega_x2))*Cross_a_x'*C_q_w_i_hat; % has error.. corrected?
    %check
    Qd (10 : 12, 1 : 3) =    (diag(var_nbomega)* (-dt^4/24*eye(3) - dt^5/120*Cross_omega_x - dt^6/720*Cross_omega_x2))* Cross_a_x*C_q_w_i_hat;
    %check
    Qd (13 : 15, 1 : 3) =    -diag(var_nba*dt^3/6)*C_q_w_i_hat;
    %Sebastian corrected
    Qd (7 : 9,   4 : 6) =    ((dt^2/2*eye(3) - dt^3/6*Cross_omega_x + dt^4/24*Cross_omega_x2) * diag(var_nomega) ...
                            + diag(var_nbomega)*(dt^4/8*eye(3) - dt^5/30*Cross_omega_x + dt^6/144*Cross_omega_x2))*Cross_a_x*C_q_w_i_hat; %has error
    %moved diagonal matrix 
    Qd (10 : 12, 4 : 6) =    (diag(var_nbomega)*(-dt^3/6*eye(3) - dt^4/24*Cross_omega_x - dt^5/120*Cross_omega_x2))*Cross_a_x* C_q_w_i_hat;
    %check
    Qd (13 : 15, 4 : 6) =    -diag(var_nba*dt^2/2)*C_q_w_i_hat;
    Qd (10 : 12, 7 : 9) =     diag(var_nbomega)*(-dt^2/2*eye(3) - dt^3/6*Cross_omega_x - dt^4/24*Cross_omega_x2);
    
    % Make it symmetric
    Qd = Qd + Qd';
    % Add diagonal entries
    Qd (1 : 3, 1 : 3) =                      diag(var_na*dt^3/3 ...
                                                + var_nba*dt^5/20) ...
                        - C_q_w_i_hat'*Cross_a_x*((dt^5/20*eye(3) + dt^7/504.*Cross_omega_x2)* diag(var_nomega) ...
                                           - diag(var_nbomega)*(dt^7/252 + dt^9/8640*Cross_omega_x2))*Cross_a_x*C_q_w_i_hat;
    %corrected 
    Qd (4 : 6, 4 : 6) =                      diag(var_na*dt...
                                                + var_nba*dt^3/3) ...
                        - C_q_w_i_hat'*Cross_a_x*((dt^3/3*eye(3) + dt^5/60*Cross_omega_x2) * diag(var_nomega) ...
                                           - diag(var_nbomega)*(dt^5/20 + dt^7/504*Cross_omega_x2))*Cross_a_x*C_q_w_i_hat;
     %corrected
    Qd (7 : 9, 7 : 9) =     diag(var_nomega*dt) ...
                          + diag(var_nbomega)*(dt^3/3*eye(3) + dt^5/60*Cross_omega_x2);
    %check                  
    Qd (10 : 12, 10 : 12) = diag(var_nbomega*dt);
    Qd (13 : 15, 13 : 15) = diag(var_nba*dt);
    Qd (16, 16) =           var_nL*dt;
       
end