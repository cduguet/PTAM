% Create Matrices out of numerical integration 

function [Q_d_num] = GetMatricesNumerical(C_q_w_i_hat, Cross_a_x, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL, delta_t)


syms tau

    %Construct submatrices of Fd
    %FIXME - WARNING -IMPORTANT: THIS IS A MODIFIED VERSION OF THE RESULTS
    %VIEWED IN STEPHAN WEISS' PAPER. THIS IS ACCORDING TO MY SOLUTION
    A = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (tau^2/2  *eye(3) - tau^3/6  * Cross_omega_x  + tau^4/24  * Cross_omega_x2);
    B = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (-tau^3/6 *eye(3) + tau^4/24 * Cross_omega_x  - tau^5/120 * Cross_omega_x2);
    C = C_q_w_i_hat'  * Cross_a_x * C_q_w_i_hat * (tau      *eye(3) - tau^2/2  * Cross_omega_x  + tau^3/6   * Cross_omega_x2);
    D = -A;
    E = eye(3) - tau*Cross_omega_x + tau^2/2 * Cross_omega_x2;
    F = -tau + tau^2/2 * Cross_omega_x - tau^3/6 * Cross_omega_x2;


    F_d =  [[  eye(3)     tau*eye(3)     A          B           -C_q_w_i_hat'*tau^2/2       zeros(3,7)];
            [ zeros(3)      eye(3)       C          D              -C_q_w_i_hat'*tau        zeros(3,7)];
            [ zeros(3)     zeros(3)      E          F                zeros(3)              zeros(3,7)];
            [ zeros(3)     zeros(3)   zeros(3)    eye(3)             zeros(3)              zeros(3,7)];
            [ zeros(3)     zeros(3)   zeros(3)   zeros(3)             eye(3)               zeros(3,7)];
            [ zeros(7,3)  zeros(7,3) zeros(7,3) zeros(7,3)          zeros(7,3)               eye(7)  ]];


    % Continuous time system noise covariance matrix
    Q_c = diag([var_na, var_nba, var_nomega, var_nbomega, var_nL]);
    % 
    G_c = [[zeros(3,13)];
           [-C_q_w_i_hat' zeros(3,10)];
           [zeros(3,6) -eye(3) zeros(3,4)];
           [zeros(3,9) eye(3) zeros(3,1)];
           [zeros(3) eye(3) zeros(3,7)];
           [zeros(1,12) 1];
           [zeros(6,13)]];


    % if it is the first iteration of the Optical Flow, consider the
    % integration time of the IMU also before the first image frame

    Q_d_num =    int(F_d*G_c*Q_c*G_c'*F_d',0,delta_t);
    Q_d_num = eval(Q_d_num);
    
end