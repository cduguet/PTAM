%% Construct H Optical Flow Velocity 

function Hv = ConstructHv(L_hat,C_q_i_c_hat,C_q_w_i_hat,Cross_v_w_i_hat,Cross_omega_x,Cross_p_i_c_hat,p_i_c_hat,v_w_i_hat)



    % Construct Measurement Matrix for velocity
        Hv = zeros(3,22);
        
    % IMPORTANT : CONSTRUCTED ACCORDING TO MY LAST SOLUTIONS FROM 16.08

    Hv(:,4:6) = L_hat * C_q_i_c_hat * C_q_w_i_hat;   %  * Delta_v_w_iC_q_w_i_hat
%     Hv(:,7:9) = L_hat * C_q_i_c_hat * C_q_w_i_hat * Cross_v_w_i_hat * C_q_w_i_hat';   %  * delta_theta_w_i
    Hv(:,16) = C_q_i_c_hat*(C_q_w_i_hat*v_w_i_hat' + Cross_omega_x*p_i_c_hat');  % * Delta_L
    Hv(:,17:19) = L_hat * C_q_i_c_hat * Cross_omega_x;                    % * Delta_p_i_c
    Hv(:,20:22) = L_hat * C_q_i_c_hat *(C_q_w_i_hat * Cross_v_w_i_hat * C_q_w_i_hat' - Cross_omega_x * Cross_p_i_c_hat* Cross_omega_x) * C_q_i_c_hat';  % * delta_theta_i_c
    
%Construct Hv from Z~  = Hv*X~

%     Hv(:,4:6) = L_hat * C_q_i_c_hat * C_q_w_i_hat;   %  * Delta_v_w_i
%     Hv(:,7:9) = -L_hat * C_q_i_c_hat * C_q_w_i_hat' * Cross_v_w_i_hat;   %  * delta_theta_w_i
%     Hv(:,16) =  C_q_i_c_hat*(C_q_w_i_hat*v_w_i_hat' + Cross_omega_x*p_i_c_hat');  % * Delta_L
%     Hv(:,17:19) = L_hat * C_q_i_c_hat * Cross_omega_x;                    % * Delta_p_i_c
%     Hv(:,20:22) = -L_hat * C_q_i_c_hat' *(C_q_w_i_hat' * Cross_v_w_i_hat - Cross_omega_x' * Cross_p_i_c_hat);  % * delta_theta_i_c

end