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

    dt =  delta_t*eye(3);
    var_na = diag(var_na);
    var_nba = diag(var_nba);
    var_nomega = diag(var_nomega);
    var_nbomega = diag(var_nbomega);
    
    
    C_q_w_i2 = C_q_w_i_hat' * C_q_w_i_hat;
    Cq__ax__Cq = C_q_w_i_hat' * Cross_a_x* C_q_w_i_hat;
    Cq__ax2__Cq = C_q_w_i_hat' * Cross_a_x^2* C_q_w_i_hat;
    wx = Cross_omega_x;
    wx2 = Cross_omega_x2;
    wx3 = -(wx(2,3)^2+wx(1,3)^2+wx(2,1)^2) * wx;
    wx4 = -(wx(2,3)^2+wx(1,3)^2+wx(2,1)^2) * wx2;
%     wx5 = (wx(2,3)^2+wx(1,3)^2+wx(2,1)^2) * wx;
    
    
    Qd (4 : 6,   1 : 3) =  (1/2)*(dt^4/4)*var_nba*C_q_w_i2 - (dt^2/2)*var_na* Cq__ax2__Cq - var_nomega* Cq__ax__Cq * (dt^4/8 - (dt^5*wx)/60 + (dt^6*wx^2)/144 + (dt^7*wx3)/ 1008 + (dt^8*wx4)/1152) * Cq__ax__Cq - var_nbomega* Cq__ax__Cq * (dt^6/72 - (dt^7*wx)/1008 - (dt^8*wx^2)/1920 + (dt^9*wx3)/ 2880 - (dt^10*wx4)/28800) * Cq__ax__Cq;
%     Qd (1 : 3,   4 : 6) =  (1/2)*(dt^4/4)*var_nba*C_q_w_i2 - (dt^2/2)*var_na* Cq__ax2__Cq - var_nomega* Cq__ax__Cq * (dt^4/8 + (dt^5*wx)/60 + (dt^6*wx^2)/144 - (dt^7*wx3)/ 1008 + (dt^8*wx4)/1152) * Cq__ax__Cq - var_nbomega* Cq__ax__Cq * (dt^6/72 + (dt^7*wx)/1008 - (dt^8*wx^2)/1920 - (dt^9*wx3)/ 2880 - (dt^10*wx4)/28800) * Cq__ax__Cq;
    
    Qd (7 : 9,   1 : 3) =  var_nomega*             (dt^3/6 - (dt^4*wx)/12 + (dt^5*wx^2)/40 + (dt^6*wx3)/144 + (dt^7*wx4)/336) * Cq__ax__Cq - var_nbomega             * (-(dt^5/30) + (dt^6*wx)/144 + (dt^7*wx^2)/5040 - (dt^8*wx3)/720 + (dt^9*wx4)/6480) * Cq__ax__Cq;
%     Qd (1 : 3,   7 : 9) = -var_nomega*Cq__ax__Cq * (dt^3/6 + (dt^4*wx)/12 + (dt^5*wx^2)/40 - (dt^6*wx3)/144 + (dt^7*wx4)/336)              + var_nbomega* Cq__ax__Cq * (-(dt^5/30) - (dt^6*wx)/144 + (dt^7*wx^2)/5040 + (dt^8*wx3)/720 + (dt^9*wx4)/6480);

    Qd (10 : 12, 1 : 3 )=  -var_nbomega            * (dt^4/24 + (dt^5*wx)/120 - (dt^6*wx^2)/720) * Cq__ax__Cq;
%     Qd (1 : 3  ,10 : 12)=   var_nbomega*Cq__ax__Cq * (dt^4/24 - (dt^5*wx)/120 - (dt^6*wx^2)/720);
    
    Qd (13 : 15, 1 : 3) =  -(1/6)*dt^3 *C_q_w_i_hat*var_nba;
%     Qd (1 : 3, 13 : 15) =  -(1/6)*dt^3             *var_nba*C_q_w_i_hat';
    
    %Sebastian corrected
    Qd (7 : 9,   4 : 6) =     var_nomega             * (dt^2/2 - (dt^3*wx)/6 + (dt^4*wx^2)/24 + (dt^5*wx3)/60 + (dt^6*wx4)/72)* Cq__ax__Cq + var_nbomega             * (dt^4/8 - (dt^5*wx)/60 + (dt^6*wx^2)/144 + (dt^7*wx3)/1008 + (dt^8*wx4)/1152) * Cq__ax__Cq;
%     Qd (4 : 6,   7 : 9) =    -var_nomega* Cq__ax__Cq * (dt^2/2 + (dt^3*wx)/6 + (dt^4*wx^2)/24 - (dt^5*wx3)/60 + (dt^6*wx4)/72)             - var_nbomega* Cq__ax__Cq * (dt^4/8 + (dt^5*wx)/60 + (dt^6*wx^2)/144 - (dt^7*wx3)/1008 + (dt^8*wx4)/1152);
    %moved diagonal matrix 
    Qd (10 : 12, 4 : 6) =    -var_nbomega            * (dt^3/6 + (dt^4*wx)/24 + (dt^5*wx^2)/120) *  Cq__ax__Cq;
%                               var_nbomega*Cq__ax__Cq * (dt^3/6 - (dt^4*wx)/24 + (dt^5*wx^2)/120);
    %check
    Qd (13 : 15, 4 : 6) =    (-(1/2))*dt^2*C_q_w_i_hat*var_nba;
%                              (-(1/2))*dt^2            *var_nba*C_q_w_i_hat';
    
    
    Qd (10 : 12, 7 : 9) =    (-(1/2))*dt^2*var_nbomega - (1/6)*dt^3*var_nbomega*wx - (1/24)*dt^4*var_nbomega*wx^2;
%                              (-(1/2))*dt^2*var_nbomega + (1/6)*dt^3*var_nbomega*wx - (1/24)*dt^4*var_nbomega*wx^2;
    % Make it symmetric
    Qd = Qd + Qd';

    % Add diagonal entries
    Qd (1 : 3, 1 : 3) =  (1/4)*(dt^5/5)*var_nba*C_q_w_i2 - (dt^3/3)*var_na*Cq__ax2__Cq - var_nomega*Cq__ax__Cq * (dt^5/20 + (dt^7*wx^2)/504 + (dt^9*wx4)/5184) * Cq__ax__Cq + var_nbomega*Cq__ax__Cq * (-(dt^7/252) + (13*dt^9*wx^2)/25920 - (dt^11*wx4)/158400) * Cq__ax__Cq;
    
    %corrected 
    Qd (4 : 6, 4 : 6) =   (dt^3/3)*var_nba*C_q_w_i2 - var_na*Cq__ax2__Cq*dt - var_nomega*Cq__ax__Cq* (dt^3/3 + (dt^5*wx^2)/60 + (dt^7*wx4)/252) * Cq__ax__Cq - var_nbomega* Cq__ax__Cq * (dt^5/20 + (dt^7*wx^2)/504 + (dt^9*wx4)/5184) * Cq__ax__Cq;
     %corrected
    Qd (7 : 9, 7 : 9) =   var_nomega*(dt + (dt^5*wx4)/20) - var_nbomega*(-(dt^3/3) - (dt^5*wx^2)/60 - (dt^7*wx4)/252);
    %check                  
    Qd (10 : 12, 10 : 12) = dt*var_nbomega; %diag(var_nbomega*dt) 
    Qd (13 : 15, 13 : 15) = dt*var_nba; %diag(var_nba*dt);
    Qd (16, 16) =           dt(1,1)*var_nL; %var_nL*dt;
    
    
    
       
       
end