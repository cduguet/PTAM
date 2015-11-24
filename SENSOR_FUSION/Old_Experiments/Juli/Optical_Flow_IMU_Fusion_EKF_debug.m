%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical Flow
% Cristian Duguet 

clc
close all
clear all

% Gravitational Constant
G =  9.80665;


% Experiment name 
experiment_date = '20120718_2021';
% experiment_date = '20120704_1913';

%% General Parameters:
%metrics in mm
sigma2IMUacc = [1 1 1]; % in (mm/s2)^2
sigma2IMUrotvel = deg2rad([5 5 5]).^2; % in rad^2


% Focal Length (mm):
f = 6; %mm 

% IMU Update Rate:
IMU_rate = 50; %Hz

% Scaling constant to be changed 
manual_scale = 36;

% Number of image instants which do not have Optical Flow
MAX_NOFLOW = 100;


%% Read files

%If the mat files are not created out of the log files, then create them
try
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
catch ME
    ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
    ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));
    
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
end

%Analyze the frequency of messages  
delta = diff(ViconTime,1);
Viconfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

delta = diff(OpticalFlowTime);
OpticalFlowfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

% FIXME: sample rate of the camera does not seem to be 15FPS 
% display(sprintf('Mean rate of flor information: %g arrays per second \n',OpticalFlowfreq);


% Coverting ground truth into mm
ViconPose(:,1:3) = ViconPose(:,1:3)*1000;

% Pixel size 
pixel_size=3.75e-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Preparation

%----------------------IMU Construction -----------------------------------

% Build an IMU sensor based on the Vicon Pose

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthetize Vicon Position to make a square: 

% make a 1x1 m square in 40 sec



ViconPose = repmat( [0 0 0 0 0 0 1],4000,1);
ViconPose(:,1) = [1:1000 1000*ones(1,1000) 999:-1:0 zeros(1,1000)]';
ViconPose(:,2) = [zeros(1,1000) 1:1000 1000*ones(1,1000) 999:-1:0 ]';


ViconTime = [1:4000]' * 1e-2;

timedifference = repmat(diff(ViconTime,1),1,3);


OpticalFlowTime = [0.0667:0.0667:40]';
OpticalFlow_x = 1280 * rand(length(OpticalFlowTime),400);
OpticalFlow_y = 960 * rand(length(OpticalFlowTime),400);
OpticalFlow_u =  5 * rand(size(OpticalFlow_x));
OpticalFlow_v =  5 * rand(size(OpticalFlow_x));

cornertime = find(OpticalFlowTime>10,1)
OpticalFlow_u(1:cornertime,:) =OpticalFlow_u(1:cornertime,:) + ones(cornertime,400);
% OpticalFlow_v(1:cornertime,:) = zeros(cornertime,400);

cornertime2 = find(OpticalFlowTime>20,1);
% OpticalFlow_u(cornertime:cornertime2,:) = zeros(cornertime2-cornertime,400);
OpticalFlow_v(cornertime:cornertime2-1,:) = OpticalFlow_v(cornertime:cornertime2-1,:) + ones(cornertime2-cornertime,400);

cornertime3 = find(OpticalFlowTime>30,1);
OpticalFlow_u(cornertime2:cornertime3-1,:) =OpticalFlow_u(cornertime2:cornertime3-1,:) -ones(cornertime2-cornertime,400);
% we are overwriting the last element in each corner, but it's ok since it
% is already passed the corner 
OpticalFlow_v(cornertime3:end,:) = OpticalFlow_v(cornertime3:end,:) -ones(length(OpticalFlowTime)-cornertime3+1,400);


matched_feat = ones(size(matched_feat))*400;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the differential quaternion 
invPose = invquat(ViconPose(1:end,4:7)); %inverse quaternion
IMUrotvel = quatmult(ViconPose(2:end,4:7),invPose(1:end-1,:));
% IMUrotvel = quat_normalize(IMUrotvel);

%FIXME: Vicon sometimes deliver bad irregular rotations. Supress these
%making it equal to their previous rotation velocity
for i=1:size(IMUrotvel,1)
    if abs(IMUrotvel(i,4) - 1) > 0.1
        IMUrotvel(i,:) = [0 0 0 1];
        ViconPose(i+1,4:7) = ViconPose(i,4:7);
        invPose(i,:) = invPose(i-1,:);
        error('hola');
    end
end

%Take the angular velocity out of the differential quaternion. 
IMUrotvel = 2*IMUrotvel(2:end,1:3)./timedifference(2:end,:);

% --- Acceleration: 
% Calculate velocity 

vel = diff(ViconPose(:,1:3)) ./timedifference; %mm/sec
IMUacc = diff(vel) ./timedifference(2:end,:) ; %mm/sec2

%%-- debug 
IMUacc_gt = IMUacc;
IMUacc_gt = quatmult([IMUacc_gt zeros(size(IMUacc_gt,1),1)],invPose(3:end,:));
IMUacc_gt = quatmult(ViconPose(3:end,4:7),IMUacc_gt);
IMUacc_gt = IMUacc_gt(:,1:3);
%%- debug


% Add Gravity
IMUacc = IMUacc + repmat([0 0 -G*1000],size(IMUacc,1),1);

% Rotate acceleration

IMUacc = quatmult([IMUacc zeros(size(IMUacc,1),1)],invPose(3:end,:));
IMUacc = quatmult(ViconPose(3:end,4:7),IMUacc);
IMUacc = IMUacc(:,1:3);


IMUTime_instants = size(IMUacc,1);

% ---- Noise: 
%FIXME; No bias considered
% Translational Noise 
% IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigma2IMUacc,IMUTime_instants,1));

% Rotational Noise
% IMUrotvel = IMUrotvel + random('Normal', zeros(IMUTime_instants,3), repmat(sigma2IMUrotvel,IMUTime_instants,1));

%-------------------------- Time adjustments ------------------------------

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

IMUTime = double((ViconTime(3:end,1)));

% the time in which both are transmiting
time_zero =  min(OpticalFlowTime(1), IMUTime(1));

OpticalFlowTime = double(OpticalFlowTime-time_zero);
IMUTime = double(IMUTime - time_zero);



for i = 1:size(matched_feat,1)
    
matched_feat(i) = min(matched_feat(i),max_feat);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Odometry and Predicion 

% ** Considerations:
% The imu data is going to be integrated according to the time steps of the
% Optical Flow (slower)
% In the firs iteration the IMU data will be integrated until its
% sample right after the new OpticalFlow sample timestamp. 
%  *IMPORTANT: While the imu is providing data in time differences (mm/s2),
%  the Optical Flow is considering velocity as just differences, not over a
%  time lapse. This means, if we want to integrate the IMU Data to the
%  Optical Flow, we have to consider  the rotational velocity as a Rotation
%  that has been made during a time slot, or a translation during that time
%  slot.

%% Counters and memory allocations
noflow_counter =0;

max_inliers = zeros(size(matched_feat));


%Sensors Information initialization
vel_IMUOF = zeros(length(OpticalFlowTime),3);
omega_IMU = zeros(length(OpticalFlowTime),3);

% % Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
% POSE_OF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison
% POSE_IMUOF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison


X_gt = repmat([ViconPose(1,1:3) 0 0 0 ViconPose(1,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1],length(OpticalFlowTime),1);
vel_gt = repmat( [0 0 0],length(OpticalFlowTime),1);
acc_gt = repmat( [0 0 0],length(OpticalFlowTime),1);

omega_m = zeros(length(OpticalFlowTime),3);
a_m = zeros(length(OpticalFlowTime),3);


%% --------------------------KALMAN FILTER---------------------------------

% we  express the first state as X0, the expected aas X_exp, the
% state itself as X 
% The error state as X_err, the prior as X_err_prior (no need to store in
% time)

%------------------------------Noises--------------------------------------

%Acceleration noise variance mm/s2
var_na = sigma2IMUacc;
% Acceleration Bias noise variance mm/s2
sigma_nba = [0.1 0.1 0.1];        var_nba = sigma_nba.^2;

% Rotational velocity noise variance
var_nomega = sigma2IMUrotvel;
%Rotational velocity bias noise variance
sigma_nbomega = [0.1 0.1 0.1];    var_nbomega = sigma_nbomega.^2;

%scale noise
var_nL = 0.2^2; % relative to unit

% Covariance Matrix of the Optical Flow measurements
R = diag([5 5 5].^2); %mm2 
%FIXME: better representation

% Noises needed by Stephan weiss q_d calculationa lgorithm 
    var_nqvw = [0 0 0];
    var_nqci = [0 0 0];
    var_npic = [0 0 0];




%-----------------------State Initialization-------------------------------

%We assume the initial position is at one meter height
p_w_i   = ViconPose(1,1:3);   %FIXME: We are taking the ground truth initial position, for comparison
v_w_i   = [0 0 0]; % Assume zero velocity 
q_w_i   = [0 0 0 1];  % Assume no rotation at the beginning 
b_omega = [0 0 0];
b_a     = [0 0 0];
L       = 1;
p_i_c   = [0 -4 -4];   %Calibration state: Distance from IMU to camera
q_i_c   = [0 0 0 1];     %calibration state: Orientation of camera w.r.t IMU frame 
 


% For the initial state, the state vector wil be 
X0 = [ p_w_i v_w_i q_w_i b_omega b_a L p_i_c q_i_c] ;  % 24 element vector 
%Expected state vector 

%State vector itself
X = repmat(X0,length(OpticalFlowTime),1);
X_exp = repmat(X0,length(OpticalFlowTime),1);

% Error state vector
Deltap_w_i = [0 0 0];
Deltav_w_i = [0 0 0];
deltaomega_w_i = [0 0 0];
Deltab_w   = [0 0 0];
Deltab_a   = [0 0 0];
DeltaL     = 0;
Deltap_i_c = [0 0 0];
deltaomega_i_c = [0 0 0];

%And the error state 
X_err0 = zeros(1,22);
X_err = repmat(X_err0,length(OpticalFlowTime),1);

% For the initial state, the error state covariance matrix will be
P_kminusone_kminusone = eye(22);




tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First iteration

% We process first the IMU information. The received Optical Flow data is 
% ignored until the first IMU data has arrived. 

OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);

delta_t = IMUTime(IMUsecondit) - IMUTime(1);


% ----------------------IMU Strapdown algorithm----------------------------
 
% The initial velocity and position information is in X0

[X_exp(1,:),omega_m(1,:),a_m(1,:)] = IMU_Strapdown(X0,IMUrotvel(1:IMUsecondit,:),IMUacc(1:IMUsecondit,:),IMUTime(1:IMUsecondit));

% -------------------------- Optical Flow----------------------------------
% Nothing to do here

% -------------------------------------------------------------------------
% ------------ Test
% testpos_OF3D = [0 0];
% testpos_OF2D = [0 0];
% testpos_IMU = [0 0];
% % 
% figure(6)
% % 
% actualpos3  = testpos_IMU;
% testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
% % 
% plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
% hold on
% -------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%%%%Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------- Expected State Values -----------------------------

p_w_i_hat       = X_exp(1,1:3);
v_w_i_hat       = X_exp(1,4:6);
q_w_i_hat       = X_exp(1,7:10);  % Expectation fof the quaternion from world to IMU 
% q_w_i_hat = quatmult(ODO_IMUrot,X(7:10)); 
% q_w_i_hat = quat_normalize(q_w_i_hat);
b_omega_hat     = X_exp(1,11:13); % Expectation of the bias of the IMU gyro measurements
b_a_hat         = X_exp(1,14:16); % Expectation of the bias of the IMU accelerometer measurements
L_hat           = X_exp(1,17);
p_i_c_hat       = X_exp(1,18:20);
q_i_c_hat       = X_exp(1,21:24);
% FIXME: b_a_hat, b_omega_hat, p_i_c_hat and q_i_c_hat are not being changed in the
% expectancy, with respect to the state K

% --------------------Construct variables out of expected values ----------

% Cross product matrices from the unbiased measurements 
Cross_a_x = VectorToCrossMat(a_m(1,:) - b_a_hat);
Cross_omega_x = VectorToCrossMat(omega_m(1,:) - b_omega_hat);
Cross_omega_x2 = Cross_omega_x*Cross_omega_x;

%Make rotation matrix out of the quaternion 
C_q_w_i_hat = QuatToRotMat(q_w_i_hat);
C_q_i_c_hat = QuatToRotMat(q_i_c_hat);

% Make cross product matrices out of v_w_i_hat and p_i_c_hat
Cross_v_w_i_hat = VectorToCrossMat(v_w_i_hat);
Cross_p_i_c_hat = VectorToCrossMat(p_i_c_hat);

%--------------------------------------------------------------------------

% ---------------------------State Matrices--------------------------------

F_d  = ConstructFd(C_q_w_i_hat,Cross_a_x,Cross_omega_x,Cross_omega_x2,delta_t);

%%------------------- Qd Calculation: Numerical -------------------------

% [F_d,Q_d_num] = GetMatricesNumerical(C_q_w_i_hat, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)

%%-------------------------------------------------------------------------

%%------------------- Qd Calculation: Analytical --------------------------

%     Qd= ConstructQdAnalytic(delta_t,C_q_w_i_hat, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)
%%-------------------------------------------------------------------------

%%-----------------------Stephan Weiss approach ----------------------------
Qd_weiss = calc_Q(delta_t,q_w_i_hat,omega_m(1,:) - b_omega_hat,a_m(1,:) - b_a_hat,...
    var_na,var_nba,...
    var_nomega,var_nbomega,...
    var_nL,...
    var_nqvw,...
    var_nqci,...
    var_npic);

Q_d = Qd_weiss(1:22,1:22);
%%-------------------------------------------------------------------------
%-----------------------PROPAGATION/PREDICTION-----------------------------

P_k_kminusone = F_d*P_kminusone_kminusone*F_d' + Q_d;

X_err_prior = X_err0 * F_d';       %Noise is of mean 0 and variance Qd --> does not affect here then  %remember the states are horizontal
%bibliography Maybeck cap 4.9


%%-------------------------UPDATE/MEASUREMENT------------------------------

% It is important to notice that we just get measurements in the velocity.
% Then, we only write the H matrix for velocity measurements 


% Construct Measurement Matrix for velocity 

H = ConstructHv(L_hat,C_q_i_c_hat,C_q_w_i_hat,Cross_v_w_i_hat,Cross_omega_x,Cross_p_i_c_hat,p_i_c_hat,v_w_i_hat);

% Compute the residual in velocity 
r = [0 0 0]';  % For this first iteration we do not have Visual measurements 

% Covariance Matrix of the Optical Flow measurements
R; % defined in initialization of Kalman Filter

%Compute the innovation 
S = H * P_k_kminusone * H' + R;

% Compute Kalman Gain 
K  = P_k_kminusone*H'/S;

% Compute correction
X_err_correction = K*r;

% Update Coviariance Matrix 
P_k_k =  (eye(22) - K*H) * P_k_kminusone * (eye(22) - K*H)' + K*R*K';

% Update Error State Vector
X_err(1,:) = X_err_prior + X_err_correction'; % again, the state is expressed horizontal

% Update State Vector 
X(OpticalFlowstart,1:6) = X_exp(1,1:6) + X_err(1,1:6);
X(OpticalFlowstart,7:10) = quatmult([0.5*X_err(1,7:9) 1], X_exp(1,7:10));
X(OpticalFlowstart,11:20) = X_exp(1,11:20) + X_err(1,10:19);
X(OpticalFlowstart,21:24) = quatmult([0.5*X_err(1,20:22) 1], X_exp(1,21:24));
X_err(OpticalFlowstart,:) = X_err(1,:);
%FIXME: we are not using the states in time instants 2:OpticalFlowTime-1

X_exp(OpticalFlowstart,:) = X_exp(1,:);
toc
% TODO: use IMU acceleration as a orientation measurement 

%%-------------------------------------------------------------------------
%% Second iteration and So on: here we include the Optical Flow

index_end = IMUsecondit;


for i = OpticalFlowstart+1: length(OpticalFlowTime)
    tic
    %Time interval
    delta_t = OpticalFlowTime(i) - OpticalFlowTime(i-1);
    
    %% ================= IMU Strapdown Algorithm ==========================
    index_start = index_end; %use the end index of last iteration
    index_end = find(IMUTime > OpticalFlowTime(i),1); % FIXME: can crash if theres no IMU sample after the last OF sample
    
    
    % The initial velocity and position information is in X(i-1,:)
    %IMPORTANT: We are including the elements index_start(actual) =
    %index_end(previous) because wee are considering the timebetween two
    %iterations. This means, we just use IMUTime(indexstart). Anyways, it
    %is done so because of symmetry
    [X_exp(i,:),omega_m(i,:),a_m(i,:)] = IMU_Strapdown(X(i-1,:),IMUrotvel(index_start:index_end,:),IMUacc(index_start:index_end,:),IMUTime(index_start:index_end));

    %% =====================OPTICAL FLOW===================================
    ama=1/36; % FIXME : get rid of this constant
    
     %% FAKE measurement
    [vel_IMUOF(i,:)] = mean(vel(index_start+2:index_end+1,:),1);
    
    %% --------------------------SCALE RECOVERY -------------------------------
    % NOT IN USE: SCALE RESOLVED IN SENSOR FUSION FILTER
    
    %      crossprod = cross(repmat(ODO_IMUrot,1,max_inliers(i)), [x';y';z'],1);
    %Here we consider that the Optical Flow is unrotated
    %      M = [[diag(u), diag(x), -repmat(vel0(1),max_inliers(i),1)];
    %           [diag(v), diag(y), -repmat(vel0(2),max_inliers(i),1)];
    %           [diag(w), diag(z), -repmat(vel0(3),max_inliers(i),1)]];
    %
    %         [~,~,V2] = svd(M);
    %         lambda_scale = V2(:,end);
    %         manual_scale = lambda_scale(end);
    
    %FIXME: Using the unifying scale right??
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %--------------------2-D tests---------------------------------------------
    %     OpticalFlow2D
    %     Assume That the camera movement is uniform in the camera. Orientation is fixed. height is fixed.
    
    %     ODO_OpticalFlow2Dpos = [-mean(double(OpticalFlow_u(i,1:max_inliers(i))))*pixelsize*POSE(3)/f  mean(double(OpticalFlow_v(i,1:max_inliers(i))))*pixelsize*POSE(3)/f];
    %     ODO_OpticalFlow2Dpos = [-mean(u)  -mean(v)];
    
    
    %   tests movements in 2D
    %     figure(1)
    %     subplot(1,3,1)
    %     title('Epipolar Search 5DoF');
    %     actualpos = testpos_OF3D;
    %     testpos_OF3D =  testpos_OF3D + ODO_OpticalFlowpos(1:2);
    %
    %     plot([actualpos(1) testpos_OF3D(1) ],[actualpos(2) testpos_OF3D(2)])
    %     hold on[zeros(7,13)]];
    
    %     % test remove outliers
    %
    %     subplot(1,3,2)
    %     title('Flow averaging 2DoF')
    %     actualpos2 = testpos_OF2D;
    %     testpos_OF2D =  testpos_OF2D + ODO_OpticalFlow2Dpos;
    %     plot([actualpos2(1) testpos_OF2D(1) ],[actualpos2(2) testpos_OF2D(2)])
    %     hold on
    %
    %     subplot(1,3,3)
    %     title('Optical Flow')
    %     quiver(full_x, full_y, full_u, full_v, 'r')
    %     hold on
    %     quiver(x , y,  u , v, 'b')
    %     xlabel('red = discarded outliers ; blue[zeros(7,13)]];
    
    %     hold off
    %
    %     subplot(2,2,4)
    %     title('Histogram')
    %     %u
    %     [all_freq,all_value] = hist(full_u,20);
    %     [inliers_freq,inliers_value] = hist(u,all_value);
    %     stem(all_value,all_freq,'r');
    %     hold on;
    %     stem(inliers_value,inliers_freq,'b');
    %     %v
    %     [all_freq,all_value] = hist(full_v,20);
    %     [inliers_freq,inliers_value] = hist(v,[zeros(7,13)]];
    
    %     stem(all_value,-all_freq,'r');
    %     stem(inliers_value,-inliers_freq,'b');
    %     hold off;
    
    
    %     title('IMU integration')
    %     actualpos3  = testpos_IMU;
    %     testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
    %     % if abs(testpos_IMU - k(j,1:2)) > 1e-5
    %     %     error('the do not coincide');
    %     % else
    %     %     fprintf('the IMU integrations coincide!! :):) \n');
    %     % end
    %     plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
    %     hold on[zeros(7,13)]];
    
    %     w = waitforbuttonpress;
    % ---------------------------------------------------------------------
        
    toc
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Error State Covariance Matrix 
    P_kminusone_kminusone = P_k_k;    

    % --------------------- Expected State Values -----------------------------
    
    p_w_i_hat       = X_exp(i,1:3);
    v_w_i_hat       = X_exp(i,4:6);
    q_w_i_hat       = X_exp(i,7:10);  % Expectation of the quaternion from world to IMU
    % q_w_i_hat = quatmult(ODO_IMUrot,X(7:10));
    % q_w_i_hat = quat_normalize(q_w_i_hat);
    b_omega_hat     = X_exp(i,11:13); % Expectation of the bias of the IMU gyro measurements
    b_a_hat         = X_exp(i,14:16); % Expectation of the bias of the IMU accelerometer measurements
    L_hat           = X_exp(i,17);
    p_i_c_hat       = X_exp(i,18:20);
    q_i_c_hat       = X_exp(i,21:24);
    % FIXME: b_a_hat, b_omega_hat, p_i_c_hat and q_i_c_hat are not being changed in the
    % expectancy, with respect to the state K
    
    
    % --------------------Construct variables out of expected values ----------
    
    
    
    % Cross product matrices from the unbiased measurements
    Cross_a_x = VectorToCrossMat(a_m(i,:) - b_a_hat);
    Cross_omega_x = VectorToCrossMat(omega_m(i,:) - b_omega_hat);
    Cross_omega_x2 = Cross_omega_x*Cross_omega_x;
    
    %Make rotation matrix out of the quaternion
    C_q_w_i_hat = QuatToRotMat(q_w_i_hat);
    C_q_i_c_hat = QuatToRotMat(q_i_c_hat);
    
    % Make cross product matrices out of v_w_i_hat and p_i_c_hat
    Cross_v_w_i_hat = VectorToCrossMat(v_w_i_hat);
    Cross_p_i_c_hat = VectorToCrossMat(p_i_c_hat);
    
    %--------------------------------------------------------------------------
    
    
    % ---------------------------State Matrices--------------------------------
    
    F_d  = ConstructFd(C_q_w_i_hat,Cross_a_x,Cross_omega_x,Cross_omega_x2,delta_t);
    
    %%------------------- Qd Calculation: Numerical -------------------------
    
%     [Q_d_num] = GetMatricesNumerical(C_q_w_i_hat, Cross_a_x, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL, delta_t);
    
    %%-------------------------------------------------------------------------
    
    %%------------------- Qd Calculation: Analytical --------------------------
    
    %     Qd= ConstructQdAnalytic(delta_t,C_q_w_i_hat, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)

    %%-------------------------------------------------------------------------
    
    %%-----------------------Stephan Weiss approach ----------------------------
    
    Qd_weiss = calc_Q(delta_t,q_w_i_hat,omega_m(i,:) - b_omega_hat,a_m(i,:) - b_a_hat,...
        var_na,var_nba,...
        var_nomega,var_nbomega,...
        var_nL,...
        var_nqvw,...
        var_nqci,...
        var_npic);
    
    Q_d = Qd_weiss(1:22,1:22);
    
    %%-------------------------------------------------------------------------
    
    %%-------------------------------------------------------------------------
    %-----------------------PROPAGATION/PREDICTION-----------------------------
    
    P_k_kminusone = F_d*P_kminusone_kminusone*F_d' + Q_d;
    
    X_err_prior = X_err(i-1,:) * F_d';      
    %Noise is of mean 0 and variance Qd --> does not affect here then  %remember the states are horizontal
    %bibliography Maybeck cap 4.9
    
    
    %%-------------------------UPDATE/MEASUREMENT------------------------------
    
    % It is important to notice that we just get measurements in the velocity.
    % Then, we only write the H matrix for velocity measurements
    
    
    H = ConstructHv(L_hat,C_q_i_c_hat,C_q_w_i_hat,Cross_v_w_i_hat,Cross_omega_x,Cross_p_i_c_hat,p_i_c_hat,v_w_i_hat);
    
    
    %It is important noting that the resulting value from the projected
    %error state vector is an expecter measurement error. This means we
    %cannot take the difference ffrom this to the actual measurement, but
    %to the difference between the actual measurement and the expected
    %measurement value. 
    
    z_exp = L_hat *C_q_i_c_hat*(C_q_w_i_hat*v_w_i_hat' + Cross_omega_x*p_i_c_hat'); 
%     z_exp =0;
    %Compute the expected mean measurement.
    
    
    % Compute the residual in velocity
    r = vel_IMUOF(i,:)' - z_exp - H * X_err_prior'; 
    
    % Covariance Matrix of the Optical Flow measurements
    R; % defined in initialization of Kalman Filter
    
    %Compute the innovation
    S = H * P_k_kminusone * H' + R;
    
    % Compute Kalman Gain
    K  = P_k_kminusone*H'/S;
    
    % Compute correction
    X_err_correction = K*r;
    
    % Update Coviariance Matrix
    P_k_k =  (eye(22) - K*H) * P_k_kminusone * (eye(22) - K*H)' + K*R*K';
    
    % Update Error State Vector
    X_err(i,:) = X_err_prior + X_err_correction'; % again, the state is expressed horizontal
    
    %FIXME: this is a debug: Ground truth
    acc_gt(i,1:3) = mean(IMUacc_gt(index_start+1:index_end,:),1);
    vel_gt(i,1:3) = mean(vel(index_start+2:index_end+1,:),1);

    X_gt(i,:) = [ViconPose(index_end+2,1:3) vel_gt(i,1:3) ViconPose(index_end+2,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1];
        
    rotvel_gt = mean(IMUrotvel(index_start+1:index_end,:),1);
    
    % Update State Vector
    X(i,1:6) = X_exp(i,1:6) + X_err(i,1:6);
    X(i,7:10) = quatmult([0.5*X_err(i,7:9) 1], X_exp(i,7:10));
    X(i,11:20) = X_exp(i,11:20) + X_err(i,10:19);
    X(i,21:24) = quatmult([0.5*X_err(i,20:22) 1], X_exp(i,21:24));
    
    toc
    % TODO: use IMU acceleration as a orientation measurement
    
end



%-------------------------Final Tests----------------------------------------
%% 3D trajectory

hPlot = figure('Name', sprintf('Position of Optical Flow'));
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
hold on
plot3(X(:,1),X(:,2),X(:,3),'r');
legend('Ground truth','EKF');
title('Optical Flow: Position');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis equal; grid on






%% Roll Pitch Yaw 
hPlot = figure('Name', sprintf('Orientation'));

[roll_EKF,pitch_EKF,yaw_EKF]=QuatToEuler(X(:,7:10));
[roll,pitch,yaw] = QuatToEuler(ViconPose(3:end,4:7));

% 
% roll_OF  =(POSE_OF(:,4));
% pitch_OF =(POSE_OF(:,5));
% yaw_OF  =(POSE_OF(:,6));
% 
% roll_IMUOF=(POSE_IMUOF(:,4));
% pitch_IMUOF=(POSE_IMUOF(:,5));
% yaw_IMUOF  =(POSE_IMUOF(:,6));
% 
% roll= (ViconPose(3:end,4));
% pitch =ViconPose(3:end,5);
% yaw = (ViconPose(3:end,4:6));



subplot(3,1,1);
title('Orientation')
plot(OpticalFlowTime(OpticalFlowstart:end),unwrap(roll_EKF(OpticalFlowstart:end)));
hold on;
plot(IMUTime,roll,'k')
xlabel('seconds'); ylabel('roll');
legend('EKF','Ground truth');
grid on

subplot(3,1,2);

plot(OpticalFlowTime(OpticalFlowstart:end),pitch_EKF(OpticalFlowstart:end));
hold on;
plot(IMUTime,pitch,'k')
xlabel('seconds'); ylabel('pitch');
legend('EKF','Ground truth');
grid on

subplot(3,1,3);

plot(OpticalFlowTime(OpticalFlowstart:end),yaw_EKF(OpticalFlowstart:end));
hold on;
plot(IMUTime,unwrap(yaw),'k')
xlabel('seconds'); ylabel('yaw');
legend('EKF','Ground truth');
grid on




%% Unrotated Linear Velocities

vel_truth =  diff(ViconPose(2:end,1:3))./timedifference(2:end,:);

hPlot = figure('Name', sprintf('Velocity of Optical Flow in world intertial frame'));
subplot(3,1,1);
title('Velocity of Optical Flow in world intertial frame');
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,1));grid on
hold on;
plot(IMUTime, vel_truth(:,1),'k');grid on
xlabel('seconds'); ylabel('velocity in x [mm/s]');
legend('EKF','Ground truth');

subplot(3,1,2);
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,2));grid on
hold on;
plot(IMUTime, vel_truth(:,2),'k');grid on
xlabel('seconds'); ylabel('velocity in y [mm/s]');
legend('EKF','Ground truth');
subplot(3,1,3);
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,3));grid on
hold on;
plot(IMUTime, vel_truth(:,3),'k');grid on
xlabel('seconds'); ylabel('velocity in z [mm/s]');
legend('EKF','Ground truth');


%% Plot scale 

hPlot = figure('Name', sprintf('Scale Evolution'));
plot(OpticalFlowTime,X(:,17))


%% Rotation velocity

hPlot = figure('Name', sprintf('Rotation Velocity'));

title('Rotation Velocity')
subplot(3,1,1);
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,1)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,1)./OFTimedifference,'r');grid on
plot(IMUTime,IMUrotvel(:,1),'k');grid on
xlabel('seconds'); ylabel('\omega_x [rad/s]');
legend('5DoF','IMU and 3DoF','Ground truth');
subplot(3,1,2);
plot(OpticalFlowTime(OpticalFlowstart:end),omega5DoF(OpticalFlowstart:end,2)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,2)./OFTimedifference,'r');grid on
plot(IMUTime,IMUrotvel(:,2),'k');grid on
xlabel('seconds'); ylabel('\omega_y [rad/s]');
legend('5DoF','IMU and 3DoF','Ground truth');
subplot(3,1,3);
plot(OpticalFlowTime(OpticalFlowstart:end),omega5DoF(OpticalFlowstart:end,3)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,3)./OFTimedifference,'r');grid on
plot(IMUTime,IMUrotvel(:,3),'k'); grid on
xlabel('seconds'); ylabel('\omega_z [rad/s]');
legend('5DoF','IMU and 3DoF','Ground truth');
grid on

subplot(3,1,1);
title('Rotation Velocity of Optical Flow')

