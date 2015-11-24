function mult = FGQGF(tau,C_q_w_i_hat,Cross_a_x,Cross_omega_x,Q_c,G_c)


A = -C_q_w_i_hat'  * Cross_a_x * (tau.^2/2  .*eye(3) - tau.^3/6  .* Cross_omega_x  + tau.^4/24  .* Cross_omega_x ^2);
B = -C_q_w_i_hat'  * Cross_a_x * (-tau.^3/6 .*eye(3) + tau.^4/24 .* Cross_omega_x  - tau.^5/120 .* Cross_omega_x ^2);
C = -C_q_w_i_hat'  * Cross_a_x * (tau      .*eye(3) - tau.^2/2  .* Cross_omega_x  + tau.^3/6   .* Cross_omega_x ^2);

D = -A;
E = eye(3) - tau.*Cross_omega_x + tau.^2/2 .* Cross_omega_x2;
F = -tau + tau.^2/2 .* Cross_omega_x - tau.^3/6 .* Cross_omega_x2;


F_d =  [[  eye(3)     tau.*eye(3)     A          B           -C_q_w_i_hat'.*tau^2/2       zeros(3,7)];
        [ zeros(3)      eye(3)       C          D              -C_q_w_i_hat'.*tau        zeros(3,7)];
        [ zeros(3)     zeros(3)      E          F                zeros(3)              zeros(3,7)];
        [ zeros(3)     zeros(3)   zeros(3)    eye(3)             zeros(3)              zeros(3,7)];
        [ zeros(3)     zeros(3)   zeros(3)   zeros(3)             eye(3)               zeros(3,7)];
        [ zeros(7,3)  zeros(7,3) zeros(7,3) zeros(7,3)          zeros(7,3)               eye(7)  ]];
  
    
    mult = F_d*G_c*Q_c*G_c'*F_d';
