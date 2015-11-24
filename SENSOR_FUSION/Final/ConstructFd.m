%% Construct Discrete Time State Transition Matrix Fd 

function F_d  = ConstructFd(C_q_w_i_hat,Cross_a_x,Cross_omega_x,Cross_omega_x2,delta_t)

%Continuous time State Transition Matrix
%     F_c = [[ zeros(3)  eye(3)           zeros(3)           zeros(3)       zeros(3)   zeros(3,7)];
%         [ zeros(3) zeros(3)   -C_q_w_i_hat'*Cross_a_x     zeros(3)   -C_q_w_i_hat'    zeros(3,7)];
%         [ zeros(3) zeros(3)       -Cross_omega_x         -eye(3)       zeros(3)    zeros(3,7)];
%         [ zeros(13,22)]];
    
    %Construct submatrices of Fd
    %FIXME - WARNING -IMPORTANT: THIS IS A MODIFIED VERSION OF THE RESULTS
    %VIEWED IN STEPHAN WEISS' PAPER. THIS IS ACCORDING TO MY SOLUTION
    A = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (delta_t^2/2  *eye(3) - delta_t^3/6  * Cross_omega_x  + delta_t^4/24  * Cross_omega_x2);
    B = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (-delta_t^3/6 *eye(3) + delta_t^4/24 * Cross_omega_x  - delta_t^5/120 * Cross_omega_x2);
    C = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (delta_t      *eye(3) - delta_t^2/2  * Cross_omega_x  + delta_t^3/6   * Cross_omega_x2);
    D = -A;
    E = eye(3) - delta_t*Cross_omega_x + delta_t^2/2 * Cross_omega_x2;
    F = -delta_t + delta_t^2/2 * Cross_omega_x - delta_t^3/6 * Cross_omega_x2;
    
%     A = -C_q_w_i_hat'  * Cross_a_x *  (delta_t^2/2  *eye(3) - delta_t^3/6  * Cross_omega_x  + delta_t^4/24  * Cross_omega_x2);
%     B = -C_q_w_i_hat'  * Cross_a_x *  (-delta_t^3/6 *eye(3) + delta_t^4/24 * Cross_omega_x  - delta_t^5/120 * Cross_omega_x2);
%     C = -C_q_w_i_hat'  * Cross_a_x *  (delta_t      *eye(3) - delta_t^2/2  * Cross_omega_x  + delta_t^3/6   * Cross_omega_x2);
%     D = -A;
%     E = eye(3) - delta_t*Cross_omega_x + delta_t^2/2 * Cross_omega_x2;
%     F = -delta_t + delta_t^2/2 * Cross_omega_x - delta_t^3/6 * Cross_omega_x2;
    
    
    
    F_d =  [[  eye(3)     delta_t*eye(3)     A          B           -C_q_w_i_hat'*delta_t^2/2       zeros(3,7)];
        [ zeros(3)      eye(3)       C          D              -C_q_w_i_hat'*delta_t        zeros(3,7)];
        [ zeros(3)     zeros(3)      E          F                zeros(3)              zeros(3,7)];
        [ zeros(3)     zeros(3)   zeros(3)    eye(3)             zeros(3)              zeros(3,7)];
        [ zeros(3)     zeros(3)   zeros(3)   zeros(3)             eye(3)               zeros(3,7)];
        [ zeros(7,3)  zeros(7,3) zeros(7,3) zeros(7,3)          zeros(7,3)               eye(7)  ]];
    
end

    % Continuous time system noise covariance matrix
%     Q_c = diag([var_na, var_nba, var_nomega, var_nbomega, var_nL]);
%     %
%     G_c = [[zeros(3,13)];
%         [-C_q_w_i_hat' zeros(3,10)];
%         [zeros(3,6) -eye(3) zeros(3,4)];
%         [zeros(3,9) eye(3) zeros(3,1)];
%         [zeros(3) eye(3) zeros(3,7)];
%         [zeros(1,12) 1];
%         [zeros(6,13)]];