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
sigma2IMUacc = (1.0e-06 *[0.1833    0.1833    0.3333]).^2; % in (mm/s2)^2
sigma2IMUrotvel = deg2rad(5e-7*[1 1 1]).^2; % in rad^2


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
timedifference = double(diff(ViconTime(:,1))) + double(diff(ViconTime(:,2)))/1e9; % sec
timedifference = repmat(timedifference,1,3);


% --- Rotation velocity: 
% The IMU is delivering the results in angular velocity


% Data Correction: We need to rotate the ViconPose by
% FIXME: Just for this experimen
%Rotation -90 in Z axis then 90 in X axis: [.5 -.5 -.5 .5]
ViconPose(:,4:7) = quatmult(repmat([.5 -.5 -.5 .5],size(ViconPose,1),1),ViconPose(:,4:7));
%FIXME: maybe its better to set the camera orientation as looking downwards
%([-1 0 0 0])

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
% IMUacc_gt = quatmult([IMUacc_gt zeros(size(IMUacc_gt,1),1)],invPose(3:end,:));
% IMUacc_gt = quatmult(ViconPose(3:end,4:7),IMUacc_gt);
% IMUacc_gt = IMUacc_gt(:,1:3);
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
IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigma2IMUacc,IMUTime_instants,1));

% Rotational Noise
IMUrotvel = IMUrotvel + random('Normal', zeros(IMUTime_instants,3), repmat(sigma2IMUrotvel,IMUTime_instants,1));

%-------------------------- Time adjustments ------------------------------

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

IMUTime = double((ViconTime(3:end,1))) + double((ViconTime(3:end,2)))/1e9 ;
OpticalFlowTime = double((OpticalFlowTime(:,1))) + double((OpticalFlowTime(:,2)))/1e9;

% the time in which both are transmiting
time_zero =  min(OpticalFlowTime(1), IMUTime(1));

OpticalFlowTime = double(OpticalFlowTime-time_zero);
IMUTime = double(IMUTime - time_zero);


% FIXME: there are still some errors in saving the files 
aux = max_feat;
max_feat = find(isnan(sum(OpticalFlow_u+OpticalFlow_v + OpticalFlow_x +OpticalFlow_y,1)),1)-1;
if isempty(max_feat), max_feat = aux;
end


OpticalFlow_u = OpticalFlow_u(:,1:max_feat);
OpticalFlow_v = OpticalFlow_v(:,1:max_feat);
OpticalFlow_x = OpticalFlow_x(:,1:max_feat);
OpticalFlow_y = OpticalFlow_y(:,1:max_feat);

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

X_gt = repmat(ViconPose(1,:),length(OpticalFlowTime),1);
vel_gt = repmat( [0 0 0],length(OpticalFlowTime),1);
acc_gt = repmat( [0 0 0],length(OpticalFlowTime),1);


omega_m = zeros(length(OpticalFlowTime),3);
a_m = zeros(length(OpticalFlowTime),3);


OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);

delta_t = IMUTime(IMUsecondit) - IMUTime(1);

X0 = [ViconPose(1,1:3) 0 0 0 ViconPose(1,4:7) 0 0 0 0 0 0]

% ----------------------IMU Strapdown algorithm----------------------------
% The initial velocity and position information is in X0

[X_exp(1,:),omega_m(1,:),a_m(1,:)] = IMU_Strapdown(X0,IMUrotvel(1:IMUsecondit,:),IMUacc(1:IMUsecondit,:),IMUTime(1:IMUsecondit));


X_exp(2:OpticalFlowstart,:) = repmat(X_exp(1,:), OpticalFlowstart-1,1);
index_end = IMUsecondit;


% fprintf('Result of the Kalman Filter: \n');
% fprintf('%d : X0_vel \t X_exp_vel \t X_err-_vel \t Measurement \t Meas_exp \t X_err_vel \t X_post_vel \t X_gt_vel \n',i);
% for ii=1:3
%     fprintf('[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
%              X0(3+ii) ,X_exp(1,3+ii), X_err_prior(3+ii),0    ,0   ,  X_err(OpticalFlowstart,3+ii)  ,    X(OpticalFlowstart,3+ii) ,  X_gt(1,3+ii) );    
% end

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
    [X_exp(i,:),omega_m(i,:),a_m(i,:)] = IMU_Strapdown([X_exp(i-1,1:16)],IMUrotvel(index_start:index_end,:),IMUacc(index_start:index_end,:),IMUTime(index_start:index_end));
 
     %FIXME: this is a debug: Ground truth
    X_gt(i,:) = ViconPose(index_end+2,:);
    acc_gt(i,1:3) = mean(IMUacc_gt(index_start+1:index_end,:),1);
    vel_gt(i,1:3) = (vel(index_end+1,:));
    
    rotvel_gt = mean(IMUrotvel(index_start+1:index_end,:),1)
end

figure('Name','acceleration comparison') 
subplot(3,1,1)
title('acceleration comparison')
plot(a_m(:,1)); hold on 
plot(acc_gt(:,1),'k'); legend('meas.','gt');
subplot(3,1,2)
plot(a_m(:,2)); hold on 
plot(acc_gt(:,2),'k'); legend('meas.','gt');
subplot(3,1,3)
plot(a_m(:,3)); hold on 
plot(acc_gt(:,3),'k'); legend('meas.','gt');


figure('Name','rotation comparison')
subplot(4,1,1)
title('rotation comparison')
plot(X_exp(:,7)); hold on 
plot(X_gt(:,4),'k'); legend('meas.','gt');
subplot(4,1,2)
plot(X_exp(:,8)); hold on 
plot(X_gt(:,5),'k'); legend('meas.','gt');
subplot(4,1,3)
plot(X_exp(:,9)); hold on 
plot(X_gt(:,6),'k'); legend('meas.','gt');
subplot(4,1,4)
plot(X_exp(:,10)); hold on 
plot(X_gt(:,7),'k'); legend('meas.','gt');



figure('Name', 'Velocity comparison')
subplot(3,1,1)
title('velocity comparison')
plot(X_exp(:,4)); hold on 
plot(vel_gt(:,1),'k'); legend('meas.','gt');
subplot(3,1,2)
plot(X_exp(:,5)); hold on 
plot(vel_gt(:,2),'k'); legend('meas.','gt');
subplot(3,1,3)
plot(X_exp(:,6)); hold on 
plot(vel_gt(:,3),'k'); legend('meas.','gt');


figure('Name','Position Comparison')
subplot(3,1,1)
title('Position comparison')
plot(X_exp(:,1)); hold on 
plot(X_gt(:,1),'k'); 
legend('meas.','gt');
subplot(3,1,2)
plot(X_exp(:,2)); hold on 
plot(X_gt(:,2),'k'); 
legend('meas.','gt');
subplot(3,1,3)
plot(X_exp(:,3)); hold on 
plot(X_gt(:,3),'k'); 
legend('meas.','gt');
