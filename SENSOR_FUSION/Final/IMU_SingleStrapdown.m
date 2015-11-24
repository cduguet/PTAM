% ----------------------IMU Strapdown algorith-----------------------------
% Calculate Odometry using the Strapdown Algorithm 
% We use the strapdown inertial navigation algorithm described in: 
% An introduction to inertial navigation - Oliver J. Woodman
% Notice that this function only affects the state elements 1:10

function X_exp = IMU_SingleStrapdown(Init_State,gyro,acc,delta_t)

% Gravitational Constant
G =  9.80665;

% steps = length(IMUTime);

X_exp = Init_State;
v_w_i_prev = Init_State(4:6)';
p_w_i_prev = Init_State(1:3)';
q_w_i_prev = quat_normalize(Init_State(7:10))'; 
C_q_w_i_prev = QuatToRotMat(q_w_i_prev);
 
b_w = Init_State(11:13)';
b_a = Init_State(14:16)';

%% Zeroth order integration

    % Rotation
    %Using  L'Hopital and the small angle assumption
    a = delta_t/2*(gyro(2,1) - b_w(1));
    b = delta_t/2*(gyro(2,2) - b_w(2));
    c = delta_t/2*(gyro(2,3) - b_w(3));
    
    q_w_i =      [[ 1  c -b a];
                  [-c  1  a b];
                  [ b -a  1 c];
                  [-a -b -c 1]] * q_w_i_prev;
            
    %NOTE: this calculation method is faster than the
    %quaternion multiplication, but less precise (numerical errors)

    q_w_i = quat_normalize(q_w_i); % We use the quaternions as row vectors instead of columns

    C_q_w_i = QuatToRotMat(q_w_i);

    %Velocity          
    v_w_i = v_w_i_prev + (C_q_w_i * (acc(2,:)' - b_a ) + [ 0 0 G]') *delta_t;
    
    %Position
    p_w_i = p_w_i_prev + v_w_i * delta_t;

    X_exp(1:3)  = p_w_i';
    X_exp(4:6)  = v_w_i';
    X_exp(7:10) = q_w_i';
    

% gyromean = mean(gyro(2,:,:),1);
% accmean = (v_w_i - Init_State(4:6)' )/delta_t;
% 
% % fprintf('duration of matrix method: %d',toc);
% 
% % tic
% % ODO_IMUrot2 = quat_normalize(Init_State(7:10))
% % 
% % for j= 2:steps
% %      delta_t_IMU = IMUTime(j) - IMUTime(j-1);
% %      ODO_IMUrot2 = quatmult([gyro(j,:).*delta_t_IMU/2 1], ODO_IMUrot2);
% % end
% % fprintf('duration of quaternion method: %d',toc);
% 
% gyromean = 2 *q_w_i(1:3)/delta_t; %FIXME: Again, small angle assumption 
% %FIXME: It does not make sense to make the small angle assumption duting
% %delta_t but not during delta_t_IMU. Why do we use then Matrix
% %multiplication? 
% 
% % Compute expectation of Rotation state (X(1,7:10) and its inverse
% X_exp(7:10) = q_w_i;
% % inv_X_exp_rot =  invquat(X_exp(7:10));
% 
% % Project accelerations onto global axis 
% % aux = quatmult([acc zeros(steps,1)], repmat(X_exp(7:10),steps,1));
% % aux = quatmult(repmat(inv_X_exp_rot,steps,1),aux);
% % ODO_IMUacc = aux(:,1:3);
% 
% % Substract gravity  component
% ODO_IMUacc(:,3) = ODO_IMUacc(:,3) + G*1000; 
% 
% % %FIXME: Not optimal at all.... just derotating back to it's orientationn in
% % %the IMU frame. 
% % accmean = quatmult([mean(ODO_IMUacc,1) 0], inv_X_exp_rot);
% % accmean = quatmult(X_exp(7:10),  accmean);
% % accmean = accmean(1:3); 
% 


% %% First order integration
% 
% a =  (gyro(1,1)  -b_w(1));
% b =  (gyro(1,2)  -b_w(2));
% c =  (gyro(1,3)  -b_w(3));
% 
% Omega_prev = [[ 1  c -b a];
%               [-c  1  a b];
%               [ b -a  1 c];
%               [-a -b -c 1]];
%               
%     % Rotation
%     
% a =  (gyro(2,1)  -b_w(1));
% b =  (gyro(2,2)  -b_w(2));
% c =  (gyro(2,3)  -b_w(3));
% 
% Omega =      [[ 1  c -b a];
%               [-c  1  a b];
%               [ b -a  1 c];
%               [-a -b -c 1]];
%               
%    
% Omega_mean = delta_t/2 * (Omega + Omega_prev)/2;
%             
% %calculate an approximated exponnential. Use minimal no. of
% %calculations
% Exponential =  eye(4);
% factorial_no =1;
% for kk = 1:5
%     factorial_no = factorial_no * kk;
%     Exponential = Exponential + Omega_mean/factorial_no;
%     Omega_mean = Omega_mean * Omega_mean;
% end
%     
% q_w_i = (Exponential + 1/48 *(Omega*Omega_prev - Omega_prev*Omega)*delta_t^2) * q_w_i_prev;
%     
% q_w_i = quat_normalize(q_w_i); % We use the quaternions as row vectors instead of columns
% 
% C_q_w_i = QuatToRotMat(q_w_i);
% 
% %Velocity         
% v_w_i = v_w_i_prev + ((C_q_w_i * (acc(2,:)' - b_a ) + C_q_w_i_prev * (acc(1,:)' -b_a ))/2 + [ 0 0 G]') *delta_t;
% 
% %Position
% p_w_i = p_w_i_prev + (v_w_i + v_w_i_prev)/2 * delta_t;
% 
% X_exp(1:3)  = p_w_i';
% X_exp(4:6)  = v_w_i';
% X_exp(7:10) = q_w_i';

% load debug.mat
% first iteration
% vel_gt = diff(ViconPose(3:11,1:3))./repmat(diff(IMUTime(1:9),1),1,3)
% acc_gt = diff(diff(ViconPose(2:11,1:3)))./repmat(diff(IMUTime(1:9),1),1,3).^2
% second iteration
% vel_gt = diff(ViconPose(11:20,1:3))./repmat(diff(IMUTime(1:10),1),1,3)
% acc_gt = diff(diff(ViconPose(10:20,1:3)))./repmat(diff(IMUTime(1:10),1),1,3).^2


%%-- debug
%acceleration without gravity component and w.r.t world frame 
% load debug.mat
% plot(diff(diff(ViconPose(2:11,1:3)))./repmat(diff(IMUTime(1:9),1),1,3).^2)
% hold on 
% plot(ODO_IMUacc(2:end,:))
end
