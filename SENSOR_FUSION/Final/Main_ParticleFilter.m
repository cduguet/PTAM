%Particle Filter attempt no. 1
%changes: 
% Units in meters instead of mm

clc
% close all
clear all

% Gravitational Constant
G =  9.80665; %m/s2

TRUE =1;
FALSE =0;


% Experiment name
% experiment_date = '20120913_1545'; % replayed experiment
% experiment_date = '20120704_1913';
% experiment_date = '20120903_stones';

experiment_date = '20120927_1857'; % Calobrated experiment - square

q_vicon_c = [0 0 -sin(pi/8) cos(pi/8)];
p_vicon_c = [-sqrt(.5)/2*3 sqrt(.5)/2 -12.5]'*1e-2; %cm->m

% q_i_c = [1 0 0 0];
% p_i_c = [.5 1.5 5.5]'*1e-2;

q_i_c = [0 0 0 1];
p_i_c = [.5 -1.5 -5.5]'*1e-2;


q_vicon_i = quatmult(invquat(q_i_c),q_vicon_c);
p_vicon_i = -QuatToRotMat(invquat(q_vicon_i))*p_i_c + p_vicon_c;



%% Read files
%If the mat files are not created from the log files, then create them

try

    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_IMU.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
    load(strcat(experiment_date,'_camera_info.txt.mat'));
    
catch ME

    ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
    ReadIMU(strcat(experiment_date,'_IMU.txt'));
    ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));
    ReadCameraInfo(strcat(experiment_date,'_camera_info.txt'));  
    
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_IMU.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
    load(strcat(experiment_date,'_camera_info.txt.mat'));
end

%Analyze the frequency of messages
delta = diff(ViconTime,1);
Viconfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

delta = diff(OpticalFlowTime);
OpticalFlowfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));
% FIXME: sample rate of the camera does not seem to be 15FPS

delta = diff(IMUTime);
IMUfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

fprintf('Mean Optical Flow rate: %g Hz\n',OpticalFlowfreq);
fprintf('Mean Vicon rate: %g Hz\n',Viconfreq);
fprintf('Mean IMU rate: %g Hz\n',IMUfreq);

% sigmaIMUacc = AngularVelocityCovariance;
% sigmaIMUrotvel = AccelerationCovariance; 

%%Observed anomaly in Imu sensor. accelerometer is rotated 180 in x- axis
%%compared to gyro

IMUrotvel(:,2:3) = -IMUrotvel(:,2:3);
IMUacc(:,1) = -IMUacc(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Preparation

%-------------------------- Time adjustments ------------------------------

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

ViconTime = double((ViconTime(:,1))) + double((ViconTime(:,2)))/1e9 ;
ViconSysTime = double((ViconSysTime(:,1))) + double((ViconSysTime(:,2)))/1e9 ;

OFTime = double((OpticalFlowTime(:,1))) + double((OpticalFlowTime(:,2)))/1e9;
OFSysTime = double((OpticalFlowSysTime(:,1))) + double((OpticalFlowSysTime(:,2)))/1e9;

IMUTime = double((IMUTime(:,1))) + double((IMUTime(:,2)))/1e9 ;
IMUSysTime = double((IMUSysTime(:,1))) + double((IMUSysTime(:,2)))/1e9 ;

CameraInfoTime = double((CameraInfoTime(:,1))) + double((CameraInfoTime(:,2)))/1e9 ;
CameraInfoSysTime = double((CameraInfoSysTime(:,1))) + double((CameraInfoSysTime(:,2)))/1e9 ;

clear('OpticalFlowTime','OpticalFlowSysTime');


%FIXME : JUST FOR THIS EXPERIMENT 
%The clock with which IMUTime was saved was much more delayed. We have to
%correct this clock bias assuming that the delay of the IMU signals were
%about the same 
OF_delay_avg = mean(OFSysTime-CameraInfoTime(4:end));
IMU_delay_avg = mean(IMUSysTime-IMUTime);
Vicon_delay_avg = mean(ViconSysTime-ViconTime);

% Assign the Camera Info Time as the Optical flow Time (time of original
% measurement)
aux = length(OFTime);
OFTime = CameraInfoTime(end-aux+1:end);


fprintf('Mean Optical Flow delay: %g s\n',OF_delay_avg);
fprintf('Mean Vicon delay: %g s\n',Vicon_delay_avg);
fprintf('Mean IMU delay: %g s\n',IMU_delay_avg);


% Expected_delay_IMU = OF_delay_avg + 100e-3;
% Expected_delay_Vicon  = OF_delay_avg + 150e-3; 

% IMUTime = IMUTime + IMU_delay_avg - Expected_delay_IMU  ; % FIXME: Adjusted Manual Delay after observations w.r.t. SIMU (Vicon)
% ViconTime = ViconTime + Vicon_delay_avg - Expected_delay_Vicon;


% % Normalize time to 0 
time_zero =  min([OFTime(1) IMUTime(1) ViconTime(1)]);

OFTime = OFTime-time_zero;
OFSysTime = OFSysTime- time_zero;
IMUTime = IMUTime-time_zero;
IMUSysTime = IMUSysTime - time_zero;
ViconTime = ViconTime-time_zero;
ViconSysTime = ViconSysTime-time_zero;


%-----------------Simulated IMU Construction ------------------------------
[SIMUacc, SIMUrotvel,SIMUTime, SIMUSysTime] = Construct_SIMU(ViconTime,ViconSysTime,ViconPose,q_vicon_i,p_vicon_i);

%-----------------Simulated OF Ouput measurements--------------------------
% How the output of the Optical Flow mean velocity between update steps
% should look like 

[SOF,SOFrotvel,SOFTime,SOFSysTime] = ConstructSOF(ViconTime,ViconSysTime,ViconPose,OFTime,OFSysTime,q_vicon_c,p_vicon_c);
%--------------------------------------------------------------------------


%   Script made to compare IMU and SIMU values. rotation, acceleration,
% %    biases unsw.
% close all
% handle1 = figure(1)
% set(handle1,'Position',[ 0 600 1500 300])
% handle2 = figure(2)
% set(handle2,'Position',[0 0 1500 300])
% 
% handle3 = figure(3)
% set(handle3,'Position',[0 0 1500 300])


% ROTATION ANALYSIS
% SIMU_rotx = cumsum(SIMUrotvel(2:end,1) .*diff(SIMUTime));
% SIMU_roty = cumsum(SIMUrotvel(2:end,2) .*diff(SIMUTime));
% SIMU_rotz = cumsum(SIMUrotvel(2:end,3) .*diff(SIMUTime));
% 
% for j = 0:1:1000
% i = 2e-3 + mod(j*1e-3,1e-2);
% IMUrotvel2 = (IMUrotvel(:,1:3) - repmat([-0.040729 +0.00272 -6e-3],length(IMUTime),1)) ;
% IMU_rotx = cumsum(IMUrotvel2(2:end,1) .*diff(IMUTime));
% IMU_roty = cumsum(IMUrotvel2(2:end,2) .*diff(IMUTime));
% IMU_rotz = cumsum(IMUrotvel2(2:end,3) .*diff(IMUTime));
% figure(handle1)
% plot(SIMUTime(2:end),SIMU_rotx, IMUTime(2:end),IMU_rotx);
% figure(handle2)
% plot(SIMUTime(2:end),SIMU_roty, IMUTime(2:end),IMU_roty);
% figure(handle3)
% plot(SIMUTime(2:end),SIMU_rotz, IMUTime(2:end),IMU_rotz);
% 
% % legend('SIMUrotx','IMUrotx','SIMUroty','IMUroty');
% % legend('SIMUrotz','IMUrotz')
% end
% % Result is: b_omega = [-0.040729 -0.00272 6e-3]


%ACCELERATION ANALYSIS
% SIMU_velx = cumsum(SIMUacc(2:end,1) .*diff(SIMUTime));
% SIMU_vely = cumsum(SIMUacc(2:end,2) .*diff(SIMUTime));
% SIMU_velz = cumsum(SIMUacc(2:end,3) .*diff(SIMUTime));
% 
% 
% for j =0:1:1
%     i  =  -2e-3 - mod(j*1e-5,1e-3);
% IMUacc2 = (IMUacc(:,1:2) - repmat([0 0],length(IMUTime),1)) ;
% IMU_velx = cumsum(IMUacc2(2:end,1) .*diff(IMUTime));
% IMU_vely = cumsum(IMUacc2(2:end,2) .*diff(IMUTime));
% IMU_velz = cumsum(IMUacc(2:end,3) .*diff(IMUTime));
% figure(handle1)
% plot(SIMUTime(2:end),SIMU_velx, IMUTime(2:end),IMU_velx);
% % plot(SIMUTime(1:end),SIMUacc(:,1), IMUTime(1:end),IMUacc(:,1));
% legend('SIMUvelx','IMUvelx');
% 
% figure(handle2)
% plot(SIMUTime(2:end),SIMU_vely, IMUTime(2:end),IMU_vely);
% % plot(SIMUTime(:),SIMUacc(:,2), IMUTime(:),IMUacc(:,2));
% legend('SIMUvely','IMUvely');
% 
% figure(handle3)
% plot(SIMUTime(2:end),SIMU_velz, IMUTime(2:end),IMU_velz);
% % plot(SIMUTime(:),SIMUacc(:,2), IMUTime(:),IMUacc(:,2));
% legend('SIMUvely','IMUvely');
% 
% legend('SIMUrotz','IMUrotz')
% end
% % Result is: b_omega = [-0.040729 -0.00272]
% 
% return
%% Script made to compare Optical Flow velocity estimation 
% 
% pixel_size = 3.75e-6;
% f = 6e-3;% Focal length
% TimeSteps_OF = size(OpticalFlow_u,1);
% OF = zeros(TimeSteps_OF,3);
% OF_omega = OF;
% 
% q_prev = ViconPose(1,4:7);
% 
% for i=2:TimeSteps_OF
%     x= OpticalFlow_x(i,:);
%     y= OpticalFlow_y(i,:);
%     u= OpticalFlow_u(i,:);
%     v= OpticalFlow_v(i,:);
%     
%     delta_t = OFTime(i)-OFTime(i-1);
%     
%     %Find orientation 
%     idx = find(ViconTime > OFTime(i),1);
%     actualtimeVicon = ViconTime(idx);
%     
%     q = ViconPose(idx,4:7);
%     
%     Rotation = QuatToRotMat(quatmult(q,invquat(q_prev)));
%     q_prev =q;
%     
%     
%     %% =====================OPTICAL FLOW===========================
%     ama=pixel_size; % FIXME : get rid of this constant
%     % Set the pixel position referenced to the center
%     x = x'-Image_width /2+0.5;
%     y = y'-Image_height/2+0.5;
%     x = x *(ama/f); y = -y *(ama/f);
%     z = ones(size(x));
%     
%     u = u'*(ama/f);
%     v = -v'*(ama/f);
%     
%     inliers_idx = HistographicRejection(u,v);
% 
%     if length(inliers_idx) < 5
%                 warning('Inliers information too poor (%d). Avoiding outlier rejection.\n',length(inliers_idx));
%                 w = zeros(size(u));
%     else
%                 x = x(inliers_idx);
%                 y = y(inliers_idx);
%                 z = z(inliers_idx);
%                 u = u(inliers_idx);
%                 v = v(inliers_idx);
%                 w = zeros(size(u));
%     end
%     
%     
%     %%UNROTATION PART
%      unrotation_OF = (eye(3)-Rotation) * [x y z]';
%             
%             %UNROTATE
%             u = u + unrotation_OF(1,:)';
%             v = v + unrotation_OF(2,:)';
%             w = w + unrotation_OF(3,:)';
%     %%%
%     
%      %% --------------and EPIPOLAR RECONSTRUCTION ------------------
%      % Calculate using the Epipolar Constraint over the Optical Flow
% %      OF(i,:) = Continuous8point1DoF(x,y,z,u,v,w)';
%        OF(i,:) = Continuous8point2DoF(x,y,z,u,v,w)'/delta_t;
% %      [aux2,aux,~] = Continuous8point5DoF(x,y,z,u,v,w);
% %      OF(i,:) = aux'; %/delta_t
% %      OF_omega(i,:) = aux2'/delta_t;
%        
% end
% L=40;
% 
% NLdistort = 2*L;
% for i=2:TimeSteps_OF%:-1:5
%     OF(i,:) = OF(i-1,:) ...
%               + sign(((OF(i,:) - OF(i-1,:))/2)).*(1 -exp(NLdistort*-abs(-OF(i,:) + OF(i-1,:))))/NLdistort;
% end
% 
% 
% for i=TimeSteps_OF:-1:5
%     OF(i,:) = [.3 .2 .4]* OF(i-2:i,:);
% %               + 0.4*sign(((-OF(i,:) + OF(i-1,:))/2)).*(1 -exp(NLdistort*abs(-OF(i,:) + OF(i-1,:))))/NLdistort;
% end
% 
% 
% 
% hPlot = figure('Name','Optical Flow Comparison');
% set(hPlot,'Position',[ 0 0 1024 768])
% set(hPlot,'PaperPositionMode', 'auto')
% % set(hPlot,'PaperUnits','normalized');
% % set(hPlot,'PaperPosition', [0 0 1 1]);
% subplot(3,1,1)
% plot(OFTime(:),OF(:,1)*L); hold on;
% plot(SOFTime(:),SOF(:,1),'k','Linewidth',2);
% xlabel('time (s)');ylabel('vel_x (m/s)');
% legend('OF','SOF'); grid on; hold off;
% subplot(3,1,2)
% plot(OFTime(:),OF(:,2)*L);hold on;
% plot(SOFTime(:),SOF(:,2),'k','Linewidth',2);
% xlabel('time (s)');ylabel('vel_y (m/s)');
% legend('OF','SOF');grid on;hold off;
% subplot(3,1,3)
% plot(OFTime(:),OF(:,3)*L);hold on;
% plot(SOFTime(:),SOF(:,3),'k','Linewidth',2);
% xlabel('time (s)');ylabel('vel_z (m/s)');
% legend('OF','SOF');grid on;hold off;
% 
% return
% 
% diffOF =(L*OF(2:end,:)-SOF);
% sprintf('Mean Square Error:\t %d', sqrt(sum(sum(diffOF.^2)))/(length(SOF)))
% saveas(hPlot,'OFcomparison','fig');
% print(hPlot,'-depsc','OFcomparison');
% save2pdf('OFcomparison',hPlot,300);


% 
% hPlot = figure('Name','Rotation Velocity');
% set(hPlot,'Position',[ 0 0 1024 768])
% set(hPlot,'PaperPositionMode','auto')
% subplot(3,1,1)
% plot(OFTime(:),OF_omega(:,1)*L); hold on;
% plot(SOFTime(:),SOFrotvel(:,1),'k','Linewidth',2);
% xlabel('time (s)');ylabel('\omega_x (rad/s)');
% legend('OF','SOF'); grid on; hold off;
% ylim([-20 20])
% subplot(3,1,2)
% plot(OFTime(:),OF_omega(:,2)*L);hold on;
% plot(SOFTime(:),SOFrotvel(:,2),'k','Linewidth',2);
% xlabel('time (s)');ylabel('\omega_y (rad/s)');
% legend('OF','SOF');grid on;hold off;
% ylim([-20 20])
% subplot(3,1,3)
% plot(OFTime(:),OF_omega(:,3)*L);hold on;
% plot(SOFTime(:),SOFrotvel(:,3),'k','Linewidth',2);
% xlabel('time (s)');ylabel('\omega_z (rad/s)');
% legend('OF','SOF');grid on;hold off;
% ylim([-20 20])
% 
% 

% diffOF =(OF_omega(2:end,:)-SOFrotvel);
% diffOF = diffOF(5:end,:)
% sprintf('Mean Square Error:\t %d', sqrt(sum(sum(diffOF.^2)))/(length(SOF)))
% saveas(hPlot,'OFcomparison_rotvel','fig');
% print(hPlot,'-depsc','OFcomparison_rotvel');
% 
% return
% Plot OOrientation
% 





%% After analysing 2012-09-03-long, found acceleration biases: 

% b_omega = [0.0036 0.0020 0];

%Rotation velocity is rotated -90 w.r.t.  louis coordinates 
% 
% theta = deg2rad(0);
% aux = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1]';
% 
% IMUrotvel = IMUrotvel * aux;
%FIXME: for this experiment


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
% Optical Flow (slower).
% In the firs iteration the IMU data will be integrated until its
% sample right after the new OpticalFlow sample timestamp.
%  *IMPORTANT: While the imu is providing data in time differences (mm/s2),
%  the Optical Flow is considering velocity as just differences, not over a
%  time lapse. This means, if we want to integrate the IMU Data to the
%  Optical Flow, we have to consider  the rotational velocity as a Rotation
%  that has been made during a time slot, or a translation during that time
%  slot.

%% Initial calibration 
% Let's assume that the camera is steady in the first 4 seconds of the
% experiments.  This calibration is still not integrated into the Filter yet,
% for real-time autocalibration.

TimeSteps_OF = length(OFTime);
TimeSteps_IMU = length(IMUTime);

k_IMU_lastcal = find(IMUTime<4,1,'last');
k_SIMU_lastcal = find(SIMUTime<4,1,'last');

% alla = zeros(k_IMU_lastcal,3);
% 
% for k= 1:k_IMU_lastcal
%     alla(k,:) = mean(IMUrotvel(1:k,:),1)
% end
% plot(alla)
b_omega = mean(IMUrotvel(1:k_IMU_lastcal,:),1);
% Here we assume that the IMU is perfectly horizontally aligned,  so that
% all the gravity component is falling in  the z accelerometer.

b_a = mean(IMUacc(1:k_IMU_lastcal,:),1) + [0 0 G];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment Setup 


% Number of image instants which do not have Optical Flow
MAX_NOFLOW = 10;

% For Optical Flow
noflow_counter =0;

%% ------------------PARTICLE FILTER INITIALIZATION -------------------------

PF = ParticleFilter;

%--------------------------Datasheet Specification-------------------------
% Datasheet specification 
% as the noise specified in ADXRS610 [  (rad/s)/sqrt(Hz) ]  * sqrt(Hz)
% sigma_omega = deg2rad([0.05 0.05 0.05]) * sqrt(IMUfreq); 

% sigma_na = ([0.6 0.6 0.9]*1e-3 *G ) * sqrt(IMUfreq); % [ in (m/s2) / sqrt (Hz)  ]  * sqrt*(Hz)
%--------------------------Acceleration Noise------------------------------

% Acceleration noise std dev m/s2 
PF.sigma_na = [1 1 2]';

% Acceleration Bias noise std dev m/s2
PF.sigma_nba = [0.5 0.5 0.5]'*1e1;

%--------------------------Angular Velocity Noise--------------------------
% 
% % Rotational velocity noise std dev
% PF.sigma_nomega = [.01 .01 0.4]*1e3';      
% 
% %Rotational velocity bias noise std dev
% PF.sigma_nbomega = [1 1 1]'*1e-5;   
% 

% Rotational velocity noise std dev
PF.sigma_nomega = [1 1 0.3]*1e1;       %  PF.var_nomega = sigma_nomega.^2;
% PF.sigma_nomega = [1 1 0.3]*1e3;   %      PF.var_nomega = sigma_nomega.^2;

%Rotational velocity bias noise variance
PF.sigma_nbomega = [0.1 0.1 0.1]*1e-2';  %  PF.var_nbomega = sigma_nbomega.^2;


%---------------------------Measurements Noise-----------------------------

% Optical Flow measurements standard deviation
PF.sigma_OF = [10 10 80]';

%---------------------------Other Parameters-------------------------------
%Number of particles 
PF.N_particles = 4000;

PF.Image_width = Image_width;
PF.Image_height = Image_height;

%position of the camera relative to the IMU
PF.p_i_c   = p_i_c;   %Calibration state: Distance from IMU to camera

% Orientation of camera w.r.t IMU frame
% PF.q_i_c   = [0 0 0.382683 0.923880]; %rotZ: 45deg    
% PF.q_i_c   = [ 0 0 0.923880 0.382683]; %rotZ:135 deg    
PF.q_i_c = q_i_c;


% FIR filter for smoothing OF 
PF.FIR_buff = zeros(3,PF.FIR_buffsize);



%%============================INITIAL STATE================================

% The state has the components: 
% X = [ p_w_i v_w_i q_w_i b_omega b_a L]

p_w_i   = ViconPose(1,1:3);   %FIXME
v_w_i   = [0 0 0]; % Assume zero velocity  
q_w_i   = quatmult(q_vicon_i,quatmult([-1 0 0 0],ViconPose(1,4:7)));  %FIXME

%remember the rotation in the last update step
PF.q_w_i_lastOF= q_w_i;

PF.L = 25;
% b_omega = [0.0033 0.00160 0];
% b_a     = [0 0 0];
% b_omega, b_ a obtained from the previous initial autocalibration %FIXME

% Initialize Particles: 
PF.X = repmat([ p_w_i v_w_i q_w_i b_omega b_a],PF.N_particles,1);

% Spread particles along b_omega, and b_a and L 

PF.X(:,11:13) = PF.X(:,11:13) + randn(PF.N_particles,3) * diag(PF.sigma_nbomega);
PF.X(:,14:16) = PF.X(:,14:16) + randn(PF.N_particles,3) * diag(PF.sigma_nba);

% Weights
PF.weights = 1/PF.N_particles * ones(PF.N_particles,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin of the iterative algorithm


% BIG ASSUMPTION: WE ASSUME HERE THAT PROCESSING TIME IS 0 . SO THIS MEANS
% WE ARE JUST CONSIDERING THE SYSTEM DELAYS OUTSIDE OF THE FILTER.
% (PF TIME, 8POINT TIME CONSIDERED ZERO)

% SECOND ASSUMPTION: MESSAGES OF EACH TOPIC ARE ARRIVING IN CONSECUTIVE
% FORM , I.E. MEASUREMENT K ALWAYS ARRIVE BEFORE MEASUREMENT K+1

TimeSteps_OF = length(OFTime);
TimeSteps_IMU = length(IMUTime);
TimeSteps_SOF =  length(SOFTime);
TimeSteps_SIMU = TimeSteps_IMU;

%-------------------------indexes------------------------------------------
k_OF = 2;
k_IMU = 2;
k_PF = 2;

%%-------------------------------------------------------------------------
% Save the Filter result in the following vectors 
% FFilter time: 
t_PF = zeros(TimeSteps_OF + TimeSteps_IMU - 2,1);
X_est = zeros(TimeSteps_OF + TimeSteps_IMU - 2,16);

OF = zeros(TimeSteps_OF-1,3);
k_upd =1;
t_OF = zeros(TimeSteps_OF-1,1);
%%-------------------------------------------------------------------------

thres_delay = 200e-3; % 100 ms threshold tolerance 
discarded = 0;
t_real =0;

total = tic % Begin clock count 


%% Only IMU + SOF

% k_SOF = 2;
% 
% while(k_SOF<=TimeSteps_SOF && k_IMU <=TimeSteps_IMU)
%     %======================================================================
%     if SOFSysTime(k_SOF) > IMUSysTime(k_IMU)
%         
%         % DISCARD MESSAGE COMING AFTER  CALCULATION
%         % If the msg corresponds to a past measurement, and the filter
%         % already calculated for that time, discard it
%         
%         if IMUTime(k_IMU) < t_PF(k_PF-1)
%             warning('discarding message!');
%             discarded = discarded+1;
%             k_IMU = k_IMU+1;
%             continue
%         end
%         
%         % WAIT FOR MESSAGE NOT TOO LATE 
%         % Wait the treshold delay time and if something "earlier" comes,
%         % then update first instead of propagate
%         if SOFSysTime(k_SOF)-IMUSysTime(k_IMU)  <  thres_delay    &&    SOFTime(k_SOF) < IMUTime(k_IMU)
%             % if the message corresponds to a past measurement, and it is
%             % already too old, discard it. This happens in case the last
%             % measurement was also an IMU measurement.
%             if SOFTime(k_SOF) <t_PF(k_PF-1)
%                 warning('discarding message!');
%                 discarded = discarded+1;
%                 k_SOF = k_SOF+1;
%                 continue
%             end
% 
%             
%             fprintf('updating...\n');
%             PF.z_SOF  = SOF(k_SOF-1,:);
%             PF = PF.Update(OpticalFlow_x(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_y(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_u(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_v(k_SOF,1:matched_feat(k_SOF)), ...
%                            SOFTime(k_SOF) - SOFTime(k_SOF-1))    ;
%             % get the posterior mean
%             X_est(k_PF,:) = PF.X_mean;
%             % save Optical Flow
%             OF(k_upd,:) = PF.z_SOF'; 
%             k_upd = k_upd+1;
%             % time of Filter equal to last processed message timestamp
%             t_PF(k_PF) = SOFTime(k_SOF);
%             
%             % if we are processing a message earlier in time, error
%             if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
%             
%             
%             k_PF = k_PF+1;
%             %time of arrival of last message (greater than t_PF)
%             t_real = max(t_real,SOFSysTime(k_SOF));
%             
%             k_SOF = k_SOF +1;
%             PF.propagated = 0;
%             continue
%         end
%         
%         
%         fprintf('propagating...\n');
%         PF = PF.Propagate(IMUrotvel(k_IMU-1:k_IMU,:), ...
%                           IMUacc(k_IMU-1:k_IMU,:), ...
%                           IMUTime(k_IMU) - IMUTime(k_IMU-1));
%         % get the posterior mean
%         X_est(k_PF,:) = PF.X_mean;
%         % time of Filter equal to last processed message timestamp
%         t_PF(k_PF) = IMUTime(k_IMU);
%         
%         % if we are processing a message earlier in time, error
%         if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
%         
%         k_PF = k_PF+1;
%         % time of arrival of last message (greater than t_PF)
%         t_real = IMUSysTime(k_IMU);
%         
%         k_IMU = k_IMU+1;
%         PF.propagated = 1;
%     %======================================================================
%     elseif IMUSysTime(k_IMU) >= SOFSysTime(k_SOF)
%         
%         % DISCARD MESSAGE COMING AFTER  CALCULATION
%         % If the msg corresponds to a past measurement, and the filter
%         % already calculated for that time, discard it
%         if SOFTime(k_SOF) <t_PF(k_PF-1)
%             warning('discarding message!');
%             discarded = discarded+1;
%             k_SOF = k_SOF+1;
%             continue
%         end
%         
%         % WAIT FOR MESSAGE NOT TOO LATE 
%         % Wait the treshold delay time and if something "earlier" comes,
%         % then update first instead of propagate
%         if IMUSysTime(k_IMU) - SOFSysTime(k_SOF) < thres_delay && SOFTime(k_SOF) > IMUTime(k_IMU)
%             % if the message corresponds to a past measurement, and it is
%             % already too old, discard it. This happens in case the last
%             % measurement was also an OF measurement.
%             if IMUTime(k_IMU) < t_PF(k_PF-1)
%                  warning('discarding message!');
%                 discarded = discarded+1;
%                 k_IMU = k_IMU+1;
%                 continue
%             end
%             
%             fprintf('propagating...\n');
%             PF = PF.Propagate(IMUrotvel(k_IMU-1:k_IMU,:), ...
%                 IMUacc(k_IMU-1:k_IMU,:), ...
%                 IMUTime(k_IMU) - IMUTime(k_IMU-1));
%             % get the posterior mean
%             X_est(k_PF,:) = PF.X_mean;
%             % time of Filter equal to last processed message timestamp
%             t_PF(k_PF) = IMUTime(k_IMU);
%             
%             % if we are processing a message earlier in time, error
%             if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
%             
%             k_PF = k_PF+1;
%             % time of arrival of last message (greater than t_PF)
%             t_real = IMUSysTime(k_IMU);
%             
%             k_IMU = k_IMU+1;
%             PF.propagated = 1;
%             continue
%         end
%         
%         % DISCARD UPDATE IN CASE THERE HAS BEEN NO PROPAGATION
%         if PF.propagated ~= 1
%             k_SOF = k_SOF+1;
%             warning('Not yet propagated at OF Timestep %d. Discarding OF.\n',k_SOF);
%             continue
%         end
%         
%         %UPDATE
%         fprintf('updating...\n');
%         PF.z_SOF  = SOF(k_SOF-1,:);
%         PF = PF.Update(OpticalFlow_x(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_y(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_u(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_v(k_SOF,1:matched_feat(k_SOF)), ...
%             SOFTime(k_SOF) - SOFTime(k_SOF-1))    ;
%         % get the posterior mean
%         X_est(k_PF,:) = PF.X_mean;
%         % save Optical Flow
%         OF(k_upd,:) = PF.z_SOF';
%         k_upd = k_upd+1;
%         % time of Filter equal to last processed message timestamp
%         t_PF(k_PF) = SOFTime(k_SOF);
%         
%         % if we are processing a message earlier in time, error
%         if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
%         
%         
%         k_PF = k_PF+1;
%         %time of arrival of last message (greater than t_PF)
%         t_real = max(t_real,SOFSysTime(k_SOF));
%         
%         k_SOF = k_SOF +1;
%         PF.propagated = 0;
%     end
%     %======================================================================
% end
% 
% 
% %%--------------------------------------------------------------------------
% 
% return


%% SIMU+ SOF

% k_SOF =2;
% k_SIMU = 2;
% 
% while(k_SOF<=TimeSteps_SOF && k_SIMU <=TimeSteps_SIMU)
%     %======================================================================
%     if SOFSysTime(k_SOF) > SIMUSysTime(k_SIMU)
%         
%         % DISCARD MESSAGE COMING AFTER  CALCULATION
%         % If the msg corresponds to a past measurement, and the filter
%         % already calculated for that time, discard it
%         
%         if SIMUTime(k_SIMU) < t_PF(k_PF-1)
%             warning('discarding message!');
%             discarded = discarded+1;
%             k_SIMU = k_SIMU+1;
%             continue
%         end
%         
%         % WAIT FOR MESSAGE NOT TOO LATE 
%         % Wait the treshold delay time and if something "earlier" comes,
%         % then update first instead of propagate
%         if SOFSysTime(k_SOF)-SIMUSysTime(k_SIMU)  <  thres_delay    &&    SOFTime(k_SOF) < SIMUTime(k_SIMU)
%             % if the message corresponds to a past measurement, and it is
%             % already too old, discard it. This happens in case the last
%             % measurement was also an SIMU measurement.
%             if SOFTime(k_SOF) <t_PF(k_PF-1)
%                 warning('discarding message!');
%                 discarded = discarded+1;
%                 k_SOF = k_SOF+1;
%                 continue
%             end
% 
%             
%             fprintf('updating...\n');
%             PF.z_SOF  = SOF(k_SOF-1,:);
%             PF = PF.Update(OpticalFlow_x(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_y(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_u(k_SOF,1:matched_feat(k_SOF)), ...
%                            OpticalFlow_v(k_SOF,1:matched_feat(k_SOF)), ...
%                            SOFTime(k_SOF) - SOFTime(k_SOF-1))    ;
%             % get the posterior mean
%             X_est(k_PF,:) = PF.X_mean;
%             % save Optical Flow
%             OF(k_upd,:) = PF.z_SOF'; 
%             k_upd = k_upd+1;
%             % time of Filter equal to last processed message timestamp
%             t_PF(k_PF) = SOFTime(k_SOF);
%             
%             % if we are processing a message earlier in time, error
%             if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
%             
%             
%             k_PF = k_PF+1;
%             %time of arrival of last message (greater than t_PF)
%             t_real = max(t_real,SOFSysTime(k_SOF));
%             
%             k_SOF = k_SOF +1;
%             PF.propagated = 0;
%             continue
%         end
%         
%         
%         fprintf('propagating...\n');
%         PF = PF.Propagate(SIMUrotvel(k_SIMU-1:k_SIMU,:), ...
%                           SIMUacc(k_SIMU-1:k_SIMU,:), ...
%                           SIMUTime(k_SIMU) - SIMUTime(k_SIMU-1));
%         % get the posterior mean
%         X_est(k_PF,:) = PF.X_mean;
%         % time of Filter equal to last processed message timestamp
%         t_PF(k_PF) = SIMUTime(k_SIMU);
%         
%         % if we are processing a message earlier in time, error
%         if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
%         
%         k_PF = k_PF+1;
%         % time of arrival of last message (greater than t_PF)
%         t_real = SIMUSysTime(k_SIMU);
%         
%         k_SIMU = k_SIMU+1;
%         PF.propagated = 1;
%     %======================================================================
%     elseif SIMUSysTime(k_SIMU) >= SOFSysTime(k_SOF)
%         
%         % DISCARD MESSAGE COMING AFTER  CALCULATION
%         % If the msg corresponds to a past measurement, and the filter
%         % already calculated for that time, discard it
%         if SOFTime(k_SOF) <t_PF(k_PF-1)
%             warning('discarding message!');
%             discarded = discarded+1;
%             k_SOF = k_SOF+1;
%             continue
%         end
%         
%         % WAIT FOR MESSAGE NOT TOO LATE 
%         % Wait the treshold delay time and if something "earlier" comes,
%         % then update first instead of propagate
%         if SIMUSysTime(k_SIMU) - SOFSysTime(k_SOF) < thres_delay && SOFTime(k_SOF) > SIMUTime(k_SIMU)
%             % if the message corresponds to a past measurement, and it is
%             % already too old, discard it. This happens in case the last
%             % measurement was also an OF measurement.
%             if SIMUTime(k_SIMU) < t_PF(k_PF-1)
%                  warning('discarding message!');
%                 discarded = discarded+1;
%                 k_SIMU = k_SIMU+1;
%                 continue
%             end
%             
%             fprintf('propagating...\n');
%             PF = PF.Propagate(SIMUrotvel(k_SIMU-1:k_SIMU,:), ...
%                 SIMUacc(k_SIMU-1:k_SIMU,:), ...
%                 SIMUTime(k_SIMU) - SIMUTime(k_SIMU-1));
%             % get the posterior mean
%             X_est(k_PF,:) = PF.X_mean;
%             % time of Filter equal to last processed message timestamp
%             t_PF(k_PF) = SIMUTime(k_SIMU);
%             
%             % if we are processing a message earlier in time, error
%             if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
%             
%             k_PF = k_PF+1;
%             % time of arrival of last message (greater than t_PF)
%             t_real = SIMUSysTime(k_SIMU);
%             
%             k_SIMU = k_SIMU+1;
%             PF.propagated = 1;
%             continue
%         end
%         
%         % DISCARD UPDATE IN CASE THERE HAS BEEN NO PROPAGATION
%         if PF.propagated ~= 1
%             k_SOF = k_SOF+1;
%             warning('Not yet propagated at OF Timestep %d. Discarding OF.\n',k_SOF);
%             continue
%         end
%         
%         %UPDATE
%         fprintf('updating...\n');
%         PF.z_SOF  = SOF(k_SOF-1,:);
%         PF = PF.Update(OpticalFlow_x(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_y(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_u(k_SOF,1:matched_feat(k_SOF)), ...
%             OpticalFlow_v(k_SOF,1:matched_feat(k_SOF)), ...
%             SOFTime(k_SOF) - SOFTime(k_SOF-1))    ;
%         % get the posterior mean
%         X_est(k_PF,:) = PF.X_mean;
%         % save Optical Flow
%         OF(k_upd,:) = PF.z_SOF';
%         k_upd = k_upd+1;
%         % time of Filter equal to last processed message timestamp
%         t_PF(k_PF) = SOFTime(k_SOF);
%         
%         % if we are processing a message earlier in time, error
%         if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
%         
%         
%         k_PF = k_PF+1;
%         %time of arrival of last message (greater than t_PF)
%         t_real = max(t_real,SOFSysTime(k_SOF));
%         
%         k_SOF = k_SOF +1;
%         PF.propagated = 0;
%     end
%     %======================================================================
% end
% 
% 
% %%--------------------------------------------------------------------------
% 
% return


%% IMU + OF

while(k_OF<=TimeSteps_OF && k_IMU <=TimeSteps_IMU)
    %======================================================================
    if OFSysTime(k_OF) > IMUSysTime(k_IMU)
        
        
        % DISCARD MESSAGE COMING AFTER  CALCULATION
        % If the msg corresponds to a past measurement, and the filter
        % already calculated for that time, discard it
        
        if IMUTime(k_IMU) < t_PF(k_PF-1)
            warning('discarding message!');
            discarded = discarded+1;
            k_IMU = k_IMU+1;
            continue
        end
        
        % WAIT FOR MESSAGE NOT TOO LATE 
        % Wait the treshold delay time and if something "earlier" comes,
        % then update first instead of propagate
        if OFSysTime(k_OF)-IMUSysTime(k_IMU)  <  thres_delay    &&    OFTime(k_OF) < IMUTime(k_IMU)
            % if the message corresponds to a past measurement, and it is
            % already too old, discard it. This happens in case the last
            % measurement was also an IMU measurement.
            if OFTime(k_OF) <t_PF(k_PF-1)
                warning('discarding message!');
                discarded = discarded+1;
                k_OF = k_OF+1;
                continue
            end

            
            fprintf('updating...\n');
            
            PF.z_SOF  = SOF(k_OF-1,:);
            PF = PF.Update(OpticalFlow_x(k_OF,1:matched_feat(k_OF)), ...
                           OpticalFlow_y(k_OF,1:matched_feat(k_OF)), ...
                           OpticalFlow_u(k_OF,1:matched_feat(k_OF)), ...
                           OpticalFlow_v(k_OF,1:matched_feat(k_OF)), ...
                           OFTime(k_OF) - OFTime(k_OF-1))    ;
            % get the posterior mean
            X_est(k_PF,:) = PF.X_mean;
            % save Optical Flow
            OF(k_upd,:) = PF.z_OF'; 
            k_upd = k_upd+1;
            % time of Filter equal to last processed message timestamp
            t_PF(k_PF) = OFTime(k_OF);
            
            % if we are processing a message earlier in time, error
            if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
            
            
            k_PF = k_PF+1;
            %time of arrival of last message (greater than t_PF)
            t_real = max(t_real,OFSysTime(k_OF));
            
            k_OF = k_OF +1;
            PF.propagated = 0;
            continue
        end
        
        
        fprintf('propagating...\n');
        PF = PF.Propagate(IMUrotvel(k_IMU-1:k_IMU,:), ...
                          IMUacc(k_IMU-1:k_IMU,:), ...
                          IMUTime(k_IMU) - IMUTime(k_IMU-1));
        % get the posterior mean
        X_est(k_PF,:) = PF.X_mean;
        % time of Filter equal to last processed message timestamp
        t_PF(k_PF) = IMUTime(k_IMU);
        
        % if we are processing a message earlier in time, error
        if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
        
        k_PF = k_PF+1;
        % time of arrival of last message (greater than t_PF)
        t_real = IMUSysTime(k_IMU);
        
        k_IMU = k_IMU+1;
        PF.propagated = 1;
    %======================================================================
    elseif IMUSysTime(k_IMU) >= OFSysTime(k_OF)
        
        % DISCARD MESSAGE COMING AFTER  CALCULATION
        % If the msg corresponds to a past measurement, and the filter
        % already calculated for that time, discard it
        if OFTime(k_OF) <t_PF(k_PF-1)
            warning('discarding message!');
            discarded = discarded+1;
            k_OF = k_OF+1;
            continue
        end
        
        % WAIT FOR MESSAGE NOT TOO LATE 
        % Wait the treshold delay time and if something "earlier" comes,
        % then update first instead of propagate
        if IMUSysTime(k_IMU) - OFSysTime(k_OF) < thres_delay && OFTime(k_OF) > IMUTime(k_IMU)
            % if the message corresponds to a past measurement, and it is
            % already too old, discard it. This happens in case the last
            % measurement was also an OF measurement.
            if IMUTime(k_IMU) < t_PF(k_PF-1)
                 warning('discarding message!');
                discarded = discarded+1;
                k_IMU = k_IMU+1;
                continue
            end
            
            fprintf('propagating...\n');
            PF = PF.Propagate(IMUrotvel(k_IMU-1:k_IMU,:), ...
                IMUacc(k_IMU-1:k_IMU,:), ...
                IMUTime(k_IMU) - IMUTime(k_IMU-1));
            % get the posterior mean
            X_est(k_PF,:) = PF.X_mean;
            % time of Filter equal to last processed message timestamp
            t_PF(k_PF) = IMUTime(k_IMU);
            
            % if we are processing a message earlier in time, error
            if(t_PF(k_PF) < t_PF(k_PF-1)),          error('Inconsistency'); end
            
            k_PF = k_PF+1;
            % time of arrival of last message (greater than t_PF)
            t_real = IMUSysTime(k_IMU);
            
            k_IMU = k_IMU+1;
            PF.propagated = 1;
            continue
        end
        
        % DISCARD UPDATE IN CASE THERE HAS BEEN NO PROPAGATION
        if PF.propagated ~= 1
            k_OF = k_OF+1;
            warning('Not yet propagated at OF Timestep %d. Discarding OF.\n',k_OF);
            continue
        end
        
        %UPDATE
        fprintf('updating...\n');
        PF.z_SOF  = SOF(k_OF-1,:);
        PF = PF.Update(OpticalFlow_x(k_OF,1:matched_feat(k_OF)), ...
            OpticalFlow_y(k_OF,1:matched_feat(k_OF)), ...
            OpticalFlow_u(k_OF,1:matched_feat(k_OF)), ...
            OpticalFlow_v(k_OF,1:matched_feat(k_OF)), ...
            OFTime(k_OF) - OFTime(k_OF-1))    ;
        % get the posterior mean
        X_est(k_PF,:) = PF.X_mean;
        % save Optical Flow
        OF(k_upd,:) = PF.z_OF';
        k_upd = k_upd+1;
        % time of Filter equal to last processed message timestamp
        t_PF(k_PF) = OFTime(k_OF);
        
        % if we are processing a message earlier in time, error
        if(t_PF(k_PF) < t_PF(k_PF-1)),      error('Inconsistency'); end
        
        
        k_PF = k_PF+1;
        %time of arrival of last message (greater than t_PF)
        t_real = max(t_real,OFSysTime(k_OF));
        
        k_OF = k_OF +1;
        PF.propagated = 0;
    end
    %======================================================================
end

totaltime = toc(total);


%-------------------------Final Tests----------------------------------------
%% 3D trajectory
 
hPlot = figure('Name', sprintf('3D Position Estimation of PF'));
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3),'k');
hold on
plot3(X_est(1:k_PF-1,1),X_est(1:k_PF-1,2),-X_est(1:k_PF-1,3) + 2*ViconPose(1,3),'b');
% plot3(OF(:,1),OF(:,2),OF(:,3),'r')
legend('Ground truth','PF');
% title('Optical Flow: Position');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
axis equal; grid on

saveas(hPlot,'PF3D','fig');
print(hPlot,'-depsc','PF3D');

%% Create a pure integrated IMU  velocity and orientation estimation

% IMUpure = repmat(ViconPose(1,:),TimeSteps_IMU,1);
IMUpure = zeros(TimeSteps_IMU,10);
IMUpure(1,1:3) = ViconPose(1,1:3);
IMUpure(1,7:10) = quatmult(q_vicon_i,ViconPose(1,4:7));

for i=2:TimeSteps_IMU
    aux = IMU_SingleStrapdown([IMUpure(i-1,:) b_omega b_a],IMUrotvel(i-1:i,:),IMUacc(i-1:i,:),IMUTime(i)-IMUTime(i-1));
    IMUpure(i,:) = aux(1:10);
end

%% Roll Pitch Yaw
hPlot = figure('Name', sprintf('Orientation'));
set(hPlot,'Position',[ 0 0 1024 768])
set(hPlot,'PaperPositionMode', 'auto')

[roll_PF,pitch_PF,yaw_PF]=QuatToEuler(X_est(1:k_PF-1,7:10));
[roll,pitch,yaw] = QuatToEuler(quatmult(repmat(q_vicon_i,size(ViconTime)),ViconPose(:,4:7)));

roll = unwrap(roll);
pitch = unwrap(pitch);
yaw = unwrap(yaw);

[roll_IMU,pitch_IMU,yaw_IMU]=QuatToEuler(IMUpure(:,7:10));
roll_IMU = unwrap(roll_IMU);


% roll_OF  =(POSE_OF(:,4));
% pitch_OF =(POSE_OF(:,5));
% yaw_OF  =(POSE_OF(:,6));
%
% roll_IMUOF=(POSE_IMUOF(:,4));
% pitch_IMUOF=(POSE_IMUOF(:,5));
% yaw_IMUOF  =(POSE_IMUOF(:,6));

% roll= (ViconPose(3:end,4));
% pitch =ViconPose(3:end,5);
% yaw = (ViconPose(3:end,4:6));


subplot(3,1,1);
title('Orientation')
plot(t_PF(1:k_PF-1),unwrap(roll_PF)-pi);
hold on;
plot(ViconTime,roll+pi,'k')
plot(IMUTime,roll_IMU+pi,'r--','Linewidth',2);
xlabel('time (s)'); ylabel('roll (rad)');
legend('PF','Ground truth','IMUfilter');
grid on

subplot(3,1,2);

plot(t_PF(1:k_PF-1),pitch_PF);
hold on;
plot(ViconTime,pitch,'k')
plot(IMUTime,pitch_IMU,'r--','Linewidth',2);
xlabel('time (s)'); ylabel('pitch (rad)');
legend('PF','Ground truth','IMUfilter');
grid on

subplot(3,1,3);

plot(t_PF(1:k_PF-1),yaw_PF-pi/4);
hold on;
plot(ViconTime,unwrap(yaw),'k')
plot(IMUTime,yaw_IMU,'r--','Linewidth',2);
xlabel('time (s)'); ylabel('yaw (rad)');
legend('PF','Ground truth','IMUfilter');
grid on

saveas(hPlot,'PFOrientation','fig');
print(hPlot,'-depsc','PFOrientation');


%% quaternion :
% hPlot = figure('Name', sprintf('OrientationQuaternion'));
% 
% subplot(4,1,1);
% title('OrientationQuaternion')
% plot(t_PF(1:k_PF-1),X(1:k_PF-1,7));
% hold on;
% plot(ViconTime,ViconPose(:,4),'k')
% % plot(IMUTime,IMUrot(:,1),'r');
% xlabel('seconds'); ylabel('rx');
% legend('PF','Ground truth','IMUfilter');
% grid on
% 
% subplot(4,1,2);
% 
% plot(t_PF(1:k_PF-1),X(1:k_PF-1,8));
% hold on;
% plot(ViconTime,ViconPose(:,5),'k')
% % plot(IMUTime,IMUrot(:,2),'r');
% xlabel('seconds'); ylabel('ry');
% legend('PF','Ground truth','IMUfilter');
% grid on
% 
% subplot(4,1,3);
% 
% plot(t_PF(1:k_PF-1),X(1:k_PF-1,9));
% hold on;
% plot(ViconTime,ViconPose(:,6),'k')
% % plot(IMUTime,IMUrot(:,3),'r');
% xlabel('seconds'); ylabel('rz');
% legend('PF','Ground truth','IMUfilter');
% grid on
% 
% subplot(4,1,4);
% 
% plot(t_PF(1:k_PF-1),X(1:k_PF-1,10));
% hold on;
% plot(ViconTime,ViconPose(:,7),'k')
% % plot(IMUTime,IMUrot(:,4),'r');
% xlabel('seconds'); ylabel('w');
% legend('PF','Ground truth','IMUfilter');
% grid on

%% Velocity

% 
% 
% vel_truth =  diff(ViconPose(:,1:3))./repmat(diff(ViconTime),1,3);
% height = interp1(ViconTime(2:end), vel_truth(:,end), t_OF);
% 
% ViconPose;
% hPlot = figure('Name', sprintf('Velocity of Optical Flow in world frame'));
% subplot(3,1,1);
% title('Velocity of Optical Flow in world intertial frame');
% plot(ViconTime(2:end), vel_truth(:,1),'k');grid on
% hold on;
% plot(t_OF,OF(:,1)*PF.X(17),'r');
% xlabel('seconds'); ylabel('vel_x (m/s)');
% legend('Ground truth','Scaled OF');
% 
% subplot(3,1,2);
% 
% plot(ViconTime(2:end), vel_truth(:,2),'k');grid on
% hold on;
% plot(t_OF,OF(:,2)*PF.X(17),'r');
% xlabel('seconds'); ylabel('vel_y (m/s)');
% legend('Ground truth','Scaled OF');
% 
% subplot(3,1,3);
% grid on
% plot(ViconTime(2:end), vel_truth(:,3),'k');grid on
% hold on;
% plot(t_OF,OF(:,3)*PF.X(17),'r');
% xlabel('seconds'); ylabel('vel_z (m/s)');
% legend('Ground truth','Scaled OF');
% 
% 
% 
% hPlot = figure('Name', sprintf('Position of PF'));
% plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
% hold on
% % plot3(X(1:k_PF-1,1),X(1:k_PF-1,2),X(1:k_PF-1,3),'r');
% plot3(OF(:,1).*height,OF(:,2).*height,OF(:,3).*height,'r')
% legend('Ground truth','PF');
% title('Optical Flow: Position');
% xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (m)');
% axis equal; grid on


%% Unrotated Linear Velocities

vel_truth =  diff(ViconPose(:,1:3))./repmat(diff(ViconTime),1,3);

C_q_vicon_i = QuatToRotMat(q_vicon_i);
vel_truth = vel_truth *C_q_vicon_i';

hPlot = figure('Name', sprintf('Velocity of PF in world frame'));
set(hPlot,'Position',[ 0 0 1024 768])
set(hPlot,'PaperPositionMode', 'auto')

subplot(3,1,1);
title('Velocity of Optical Flow in world intertial frame');
plot(t_PF(1:k_PF-1),X_est(1:k_PF-1,4));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,1),'k');grid on
plot(IMUTime(2:end),IMUpure(2:end,4),'r');
xlabel('seconds'); ylabel('vel_x (m/s)');
legend('PF','Ground truth','Integrated IMU');
ylim([-1.5 1.5])

subplot(3,1,2);
plot(t_PF(1:k_PF-1),X_est(1:k_PF-1,5));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,2),'k');grid on
plot(IMUTime(2:end),IMUpure(2:end,5),'r');
xlabel('seconds'); ylabel('vel_y (m/s)');
legend('PF','Ground truth','Integrated IMU');
ylim([-1.5 1.5])

subplot(3,1,3);
plot(t_PF(1:k_PF-1),-X_est(1:k_PF-1,6));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,3),'k');grid on
plot(IMUTime(2:end),-IMUpure(2:end,6),'r');
xlabel('seconds'); ylabel('vel_z (m/s)');
legend('PF','Ground truth','Integrated IMU');
ylim([-1 1])

saveas(hPlot,'PFvelocity','fig');
print(hPlot,'-depsc','PFvelocity');

% 
%% Plot scale

hPlot = figure('Name', sprintf('Scale'));
set(hPlot,'Position',[ 0 0 300 300])
set(hPlot,'PaperPositionMode', 'auto')
plot(t_PF(1:k_PF-1),X_est(1:k_PF-1,17)); grid on;
xlabel('time (s)');% ylabel('magnitude')

saveas(hPlot,'PFscale','fig');
print(hPlot,'-depsc','PFscale');

% 
%% Plot rotation bias

hPlot = figure('Name', sprintf('Rotation velocity bias'));
set(hPlot,'Position',[ 0 0 300 300])
set(hPlot,'PaperPositionMode', 'auto')
% title('Rotation Velocity Bias');
plot(t_PF(1:k_PF-1),X_est(1:k_PF-1,11:13)); grid on;
legend('b_{\omega x}','b_{\omega y}','b_{\omega z}');
xlabel('time (s)');ylabel('b_\omega');

saveas(hPlot,'PFbomega','fig');
print(hPlot,'-depsc','PFbomega');

%% Plot acceleration bias
hPlot = figure('Name', sprintf('Acceleration bias'));
set(hPlot,'Position',[ 0 0 300 300])
set(hPlot,'PaperPositionMode', 'auto')
% title('Acceleration Bias');
plot(t_PF(1:k_PF-1),X_est(1:k_PF-1,14:16)); grid on
axis('image')
legend('b_{a x}','b_{a y}','b_{a z}');
xlabel('time (s)');ylabel('b_a');

saveas(hPlot,'PFba','fig');
print(hPlot,'-depsc','PFba');



            



% h1=plot3(X_est(:,1),X_est(:,2),X_est(:,3),'g--o'); grid on; hold on;
% h2=plot3(POSTERIOR(:,1),POSTERIOR(:,2),POSTERIOR(:,3),'-k*');
% h3=plot3(L1(1),L1(2),L1(3),'rs');
% plot3(L2(1),L2(2),L2(3),'rs');
% plot3(L3(1),L3(2),L3(3),'rs');
% plot3(L4(1),L4(2),L4(3),'rs');
% %particles
% % h4=plot3(X_P(:,1),X_P(:,2),X_P(:,3),'b.');
% legend([h1 h2 h3],'True Path','Estimated Path','Lanscapes')
% 
% 
% figure(2)
% plot(error); title('Evolution of the error on time'); grid on;
% xlabel('time (steps)'); ylabel('Error (m)')
% 
% figure(3)
% plot(Neff); title('Evolution of effective number of particles'); grid on; 
% % hold on;  plot(repmat(N_particles,1,K));
% xlabel('time (steps)'); ylabel('Neff')
% 
% %create annotations to plot
% % sorted_w = sort(w_k,'descend');
% % for j=1:N_particles
% %     annotation('textbox',[ 0.5 0.6 0.1 0.05],...
% %         'String', sprintf('%g',find(sorted_w == w_k(j))),...
% %         'FontSize',8);
% % end

