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
sigmaIMUacc = [1 1 1]; % in (mm/s2)^2
sigmaIMUrotvel = deg2rad([5 5 5]); % in rad^2


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

% Take the differential quaternion 
invPose = invquat(ViconPose(1:end-1,4:7)); %inverse quaternion
IMUrotvel = quatmult(ViconPose(2:end,4:7),invPose);
% IMUrotvel = quat_normalize(IMUrotvel);

%FIXME: Vicon sometimes deliver bad irregular rotations. Supress these
%making it equal to their previous rotation velocity
for i=1:size(IMUrotvel,1)
    if abs(IMUrotvel(i,4) - 1) > 0.1
        IMUrotvel(i,:) = [0 0 0 1];
        ViconPose(i+1,4:7) = ViconPose(i,4:7);
        invPose(i,:) = invPose(i-1,:);
    end
end
%Take the angular velocity out of the differential quaternion. 
IMUrotvel = 2*IMUrotvel(2:end,1:3)./timedifference(2:end,:);


% --- Acceleration: 
% Calculate velocity 

vel = diff(ViconPose(:,1:3)) ./timedifference; %mm/sec
IMUacc = diff(vel) ./timedifference(2:end,:) ; %mm/sec2

% Add Gravity
IMUacc = IMUacc + repmat([0 0 G*1000],size(IMUacc,1),1);

% Rotate acceleration

IMUacc = quatmult([IMUacc zeros(size(IMUacc,1),1)],invPose(1:end-1,:));
IMUacc = quatmult(ViconPose(2:end-1,4:7),IMUacc);
IMUacc = IMUacc(:,1:3);

IMUTime_instants = size(IMUacc,1);


% ---- Noise: 
%FIXME; No bias considered
% Translational Noise 
IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigmaIMUacc,IMUTime_instants,1));

% Rotational Noise
% IMUrotvel = IMUrotvel + random('Normal', zeros(IMUTime_instants,3), repmat(sigmaIMUrotvel,IMUTime_instants,1));

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
if isempty(max_feat) max_feat = aux;
end


OpticalFlow_u = OpticalFlow_u(:,1:max_feat);
OpticalFlow_v = OpticalFlow_v(:,1:max_feat);
OpticalFlow_x = OpticalFlow_x(:,1:max_feat);
OpticalFlow_y = OpticalFlow_y(:,1:max_feat);

for i = 1:size(matched_feat,1)
    
matched_feat(i) = min(matched_feat(i),max_feat);

end

% ---------- Test ---------------------------------------------------------
% figure(1)
% plot(ViconPose(:,1),ViconPose(:,2))

% Test: to test how it would be different if we know the initial
% velocity
% h = cumsum((cumsum(IMUacc .*timedifference(2:end,:),1)+repmat(vel(1,:),size(IMUacc,1),1))  .*timedifference(2:end,:), 1);
% figure(2)
% plot(h(:,1),h(:,2))
% 
% kk = cumsum(IMUacc.*timedifference(2:end,:) ,1);
% k = cumsum( kk.*timedifference(2:end,:) ,1) ;
% figure(3)
% plot(k(:,1),k(:,2))
% 
% 
% fprintf('Debugging: \t Step \t Time Difference  \t Acc \t Vel \n' );
% for i=2:23
% fprintf('\t %u | %d | %d %d | %d %d | %d %d | \n\n ',i,timedifference(i+1,1),IMUacc(i,1:2), kk(i,1:2), k(i,1:2));
% end

% Here ends the data preparation and begins the Odometry algorithm
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
vel3DoF = zeros(length(OpticalFlowTime),3);
vel5DoF = zeros(length(OpticalFlowTime),3);
omega5DoF = zeros(length(OpticalFlowTime),3);
omega_IMU = zeros(length(OpticalFlowTime),3);

% Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
POSE_OF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison
POSE_IMUOF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First iteration

% We process first the IMU information. The received Optical Flow data is 
% ignored until the first IMU data has arrived. 

OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);


% --------------------------    IMU     -----------------------------------
 
% we don't use IMUacc(0,:)

% Initial velocity: Very important!!! Otherwise the position error will be
% propagating
% ODO_IMUvel = vel(2,:);  %[0 0 0];
% ODO_IMUvel = IMUacc(1,:) .* double(ViconTime(3) - ViconTime(2)) * 1e-9;
ODO_IMUvel = [0 0 0];
% FIXME: still not including the initial zero velocity assumption to
% estimate the gyro bias

% start at odometry movement (0,0,0) and rotation (0,0,0,1) (neutral quaternion)
% ODO_IMUpos = ODO_IMUvel.* double(ViconTime(3) - ViconTime(2)) * 1e-9
ODO_IMUpos = [0 0 0];
ODO_IMUrot = [ 0 0 0 1]; % no rotation

% --- debug:
% fprintf('Debugging: \t Step \t Time Difference  \t Acc \t Vel \n' );
Thetarot  = eye(4);
for j = 2: IMUsecondit
    
    % Rotation
    %Using  L'Hopital and the small angle assumption
    delta_t_IMU = IMUTime(j) - IMUTime(j-1);
    
    Thetarot = [[            1                delta_t_IMU/2*IMUrotvel(j,3) -delta_t_IMU/2*IMUrotvel(j,2) delta_t_IMU/2*IMUrotvel(j,1)];
                [-delta_t_IMU/2*IMUrotvel(j,3)                1             delta_t_IMU/2*IMUrotvel(j,1) delta_t_IMU/2*IMUrotvel(j,2)];
                [ delta_t_IMU/2*IMUrotvel(j,2) -delta_t_IMU/2*IMUrotvel(j,1)           1                 delta_t_IMU/2*IMUrotvel(j,3)];
                [-delta_t_IMU/2*IMUrotvel(j,1) -delta_t_IMU/2*IMUrotvel(j,2) -delta_t_IMU/2*IMUrotvel(j,3)           1               ]] * Thetarot;
      
    
    ODO_IMUrot = quatmult([0.5*IMUrotvel(j,:)*delta_t_IMU 1],ODO_IMUrot);
    % Translation : NOTICE IT IS NOT YET UNROTATED!!!
    
    ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU;
    ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
end

% ODO_IMUrot = quat_normalize(ODO_IMUrot*Thetarot'); % We use the quaternions as row vectors instead of columns

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



%% Second iteration and So on: here we include the Optical Flow


index_end = IMUsecondit;

for i = OpticalFlowstart+1: length(OpticalFlowTime)
    %Time interval
    delta_t = OpticalFlowTime(i) - OpticalFlowTime(i-1);
    
    %% ========================== IMU =====================================   
    index_start = index_end; %use the end index of last iteration
    index_end = find(IMUTime > OpticalFlowTime(i),1); % FIXME: can crash if theres no IMU sample after the last OF sample 
    
    ODO_IMUpos = [0 0 0];
    ODO_IMUrot = [ 0 0 0 1]; % no rotation 
    
    for j = index_start+1: index_end    
        % Translation: These translations are recovered w.r.t the camera
        % frame
        delta_t_IMU = IMUTime(j) - IMUTime(j-1);
        ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU;
        ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
        
        
        % Rotation
        %Using  L'Hopital and the small angle assumption
        Thetarot = [[            1                delta_t_IMU/2*IMUrotvel(j,3) -delta_t_IMU/2*IMUrotvel(j,2) delta_t_IMU/2*IMUrotvel(j,1)];
                    [-delta_t_IMU/2*IMUrotvel(j,3)                1             delta_t_IMU/2*IMUrotvel(j,1) delta_t_IMU/2*IMUrotvel(j,2)];
                    [ delta_t_IMU/2*IMUrotvel(j,2) -delta_t_IMU/2*IMUrotvel(j,1)           1                 delta_t_IMU/2*IMUrotvel(j,3)];
                    [-delta_t_IMU/2*IMUrotvel(j,1) -delta_t_IMU/2*IMUrotvel(j,2) -delta_t_IMU/2*IMUrotvel(j,3)           1               ]] * Thetarot;
                
        ODO_IMUrot = quatmult([0.5*IMUrotvel(j,:)*delta_t_IMU 1],ODO_IMUrot);
    end
    
%     ODO_IMUrot = quat_normalize(ODO_IMUrot*Thetarot'); % We use the quaternions as row vectors instead of columns
    
    
    %% =====================OPTICAL FLOW===================================
    ama=1/36;
    
    % Set the pixel position referenced to the center
    x = (OpticalFlow_x(i,1:matched_feat(i)))'-Image_width/2+0.5;
    y = (OpticalFlow_y(i,1:matched_feat(i)))'-Image_height/2+0.5;
    x = x *(ama/f); y = -y *(ama/f);
    z = ones(size(x));
    
    u = (OpticalFlow_u(i,1:matched_feat(i)))'*(ama/f);
    v = -(OpticalFlow_v(i,1:matched_feat(i)))'*(ama/f);      
    
    %FIXME now they are just being divided by f, but not multiplied by
    %pixelsize
    
    
    %%In case of no flow detected
    if length(u) < 5
        if  noflow_counter > MAX_NOFLOW
            error('Optical Flow information too poor. Need more features.\n');
            return
        end
        warning('No flow detected in this image frame.\n');
        noflow_counter  = noflow_counter+1;
        %keep constant movement
        POSE_OF(i,:) = [ POSE_OF(i-1,1:3)+vel3DoF(i,:)'*manual_scale  quatmult([0.5*omega 1],POSE_OF(i-1,:))];
        continue
    end
    
       
    %% ---------------------RANSAC-----------------------------------------
    
    
    
    % First method: Calculating the Homogrpahy Matrix from a random subset
    % of sample_num points of all the matched features in the optical flow.
    % After max_Iter iterations, return the inliers of the most accepted
%     % model (The model which can include the most inliers).
%     max_Iter = 1000;
%     sample_num = 5;
%     RANSAC_Homography_treshold = 1;
%     
%     inliers_idx = RANSAC_Homography(x,y,u,v,max_Iter,sample_num,RANSAC_Homography_treshold);

    % Second method: Choosing the Consensus of the random samples in a
    % linear sense (without Homography model) 
% 
    max_Iter = 1000;
    sample_num = 50;    
    RANSAC_Linear_treshold = 1;

    inliers_idx = RANSAC_Linear(u,v,max_Iter,sample_num,RANSAC_Linear_treshold);

     if length(inliers_idx) < 5
         warning('Inliers information too poor. Avoiding RANSAC selection.\n');
         w = zeros(size(u));
     else
         full_x=x;
         full_y=y;
         full_u=u;
         full_v=v;
         x = x(inliers_idx);
         y = y(inliers_idx);
         z = z(inliers_idx);
         u = u(inliers_idx);
         v = v(inliers_idx);
         w = zeros(size(u));
         
         max_inliers(i) = length(inliers_idx);
     end


    %%----------------RANSAC_END-------------------------------------------
    


    %% ---------------EPIPOLAR RECONSTRUCTION- 5DoF-------------------------
    
    %     Case 5DoF
    [omega5DoF(i,:),vel5DoF(i,:),correl5DoF] = Continuous8point5DoF(x,y,z,u,v,w);
    
    % Case 3DoF w=(0,0,0) FIRST UNROTATE!!!
    % As we get the rotation information from the IMU, we restrict the
    % solution of the OF epipolar search to the search of the translational
    % velocity
    
    %% ------------ SENSOR COLLIGATION and EPIPOLAR RECONSTRUCTION 2DoF----
    % Unrotate Optical Flow for use in 3DoF Epipolar Algorithm
    % The rotation velocity of the IMU has been integrated over time. The
    % result of it is the rotation quaternion. As the Optical Flow here is
    % a distance vector rather than a velocity vector, we use the
    % integrated rotation of the IMU to unrotate the Optical Flow
    %
%     quatmult( [x y z 0]),ODO_IMUrot);
    % Assuming smart angle:
    
    omega_IMU(i,:) = 2* ODO_IMUrot(1:3);
    unrotation = cross(repmat(omega_IMU(i,:),max_inliers(i),1),[x y z],2);
    
    %UNROTATE
    u = u + unrotation(:,1);
    v = v + unrotation(:,2);
    w = w + unrotation(:,3);

    
    [vel3DoF(i,:)] = Continuous8point3DoF(x,y,z,u,v,w);
    
    %  Case 2DoF w=(0,0,0), v3 =0;
    %  A = [u3.*y-v.*z     u.*z-u3.*x];
    
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
    
    %% Assign Odometry
    
    % Recover angular motion and velocity
    %Convert the angular velocity
    
    %Scale recovery still does not work.... will try in EKF
    
    % Integrate Rotation Odometry to the previous state
    
    
    POSE_OF(i,4:7) = quatmult([0.5*omega5DoF(i,:) 1],POSE_OF(i-1,4:7));
    POSE_IMUOF(i,4:7) = quatmult(ODO_IMUrot, POSE_IMUOF(i-1,4:7));
    
    
    %Derotating camera measurements
    invrot = invquat(POSE_OF(i-1,4:7));
    aux = quatmult([ vel5DoF(i,:)*manual_scale 0], POSE_OF(i-1,4:7));
    aux = quatmult(invrot,aux);
    POSE_OF(i,1:3) =      POSE_OF(i-1,1:3)   +  aux(1:3);
    
    
    %IMU and 3DoF Optical Flow 
    invrot = invquat(POSE_IMUOF(i-1,4:7));
    aux = quatmult([ vel3DoF(i,:)*manual_scale 0], POSE_IMUOF(i-1,4:7));
    aux = quatmult(invrot,aux);
    POSE_IMUOF(i,1:3) =      POSE_IMUOF(i-1,1:3)   +  aux(1:3);
    
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
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %--------------------2-D tests---------------------------------------------
   
    %   tests movements in 2D
    %     figure(1)
    %     subplot(1,3,1)
    %     title('Epipolar Search 5DoF');
    %     actualpos = testpos_OF3D;
    %     testpos_OF3D =  testpos_OF3D + ODO_OpticalFlowpos(1:2);
    %
    %     plot([actualpos(1) testpos_OF3D(1) ],[actualpos(2) testpos_OF3D(2)])
    %     hold on
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
    %     xlabel('red = discarded outliers ; blue = inliers');
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
    %     [inliers_freq,inliers_value] = hist(v,all_value);
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
    %     hold on
    %     w = waitforbuttonpress;
    % ---------------------------------------------------------------------
    
    
    
    
    
    
    
end


%-------------------------Final Tests----------------------------------------
%% 3D trajectory
hPlot = figure('Name', sprintf('Position of Optical Flow'));
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
hold on
plot3(POSE_OF(:,1),POSE_OF(:,2),POSE_OF(:,3),'r');
plot3(POSE_IMUOF(:,1),POSE_IMUOF(:,2),POSE_IMUOF(:,3),'k');
legend('Ground truth','Optical Flow 5DoF','Optical Flow 3DoF');
title('Optical Flow: Position');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis equal; grid on

%% Velocities in time: Rotated
hPlot = figure('Name', sprintf('Velocity of Optical Flow'));

OFTimedifference = diff(OpticalFlowTime(OpticalFlowstart-1:end));

subplot(3,1,1);
title('Velocity of Optical Flow')
plot(OpticalFlowTime(OpticalFlowstart:end),vel5DoF(OpticalFlowstart:end,1)*manual_scale./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),vel3DoF(OpticalFlowstart:end,1)*manual_scale./OFTimedifference,'r');grid on
xlabel('seconds'); ylabel('velocity in x [mm/s]');
legend('5DoF','IMU and 3DoF');

subplot(3,1,2);
plot(OpticalFlowTime(OpticalFlowstart:end),vel5DoF(OpticalFlowstart:end,2)*manual_scale./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),vel3DoF(OpticalFlowstart:end,2)*manual_scale./OFTimedifference,'r');grid on
xlabel('seconds'); ylabel('velocity in y [mm/s]');
legend('5DoF','IMU and 3DoF');
subplot(3,1,3);
plot(OpticalFlowTime(OpticalFlowstart:end),vel5DoF(OpticalFlowstart:end,3)*manual_scale./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),vel3DoF(OpticalFlowstart:end,3)*manual_scale./OFTimedifference,'r');grid on
xlabel('seconds'); ylabel('velocity in z [mm/s]');
legend('5DoF','IMU and 3DoF');


%% Unrotated Linear Velocities
vel_truth =  diff(ViconPose(2:end,1:3))./timedifference(2:end,:);
rotatedvel_OF = (POSE_OF(OpticalFlowstart:end,:)-POSE_OF(OpticalFlowstart-1:end-1,:));
rotatedvel_IMUOF = (POSE_IMUOF(OpticalFlowstart:end,:)-POSE_IMUOF(OpticalFlowstart-1:end-1,:));

hPlot = figure('Name', sprintf('Velocity of Optical Flow in world intertial frame'));
subplot(3,1,1);
title('Velocity of Optical Flow in world intertial frame');
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_OF(:,1)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_IMUOF(:,1)./OFTimedifference,'r');grid on
plot(IMUTime, vel_truth(:,1),'k');grid on
xlabel('seconds'); ylabel('velocity in x [mm/s]');
legend('5DoF','IMU and 3DoF','Ground truth');

subplot(3,1,2);
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_OF(:,2)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_IMUOF(:,2)./OFTimedifference,'r');grid on
plot(IMUTime, vel_truth(:,2),'k');grid on
xlabel('seconds'); ylabel('velocity in y [mm/s]');
legend('5DoF','IMU and 3DoF','Ground truth');
subplot(3,1,3);
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_OF(:,3)./OFTimedifference);grid on
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),rotatedvel_IMUOF(:,3)./OFTimedifference,'r');grid on
plot(IMUTime, vel_truth(:,3),'k');grid on
xlabel('seconds'); ylabel('velocity in z [mm/s]');
legend('5DoF','IMU and 3DoF','Ground truth');



%% Rotation velocity

hPlot = figure('Name', sprintf('Rotation Velocity of Optical Flow'));

title('Rotation Velocity of Optical Flow')
subplot(3,1,1);
plot(OpticalFlowTime(OpticalFlowstart:end),omega5DoF(OpticalFlowstart:end,1)./OFTimedifference);grid on
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
%% Roll Pitch Yaw 
hPlot = figure('Name', sprintf('Orientation of Optical Flow'));

[roll_OF,pitch_OF,yaw_OF]  =QuatToEuler(POSE_OF(:,4:7));
[roll_IMUOF,pitch_IMUOF,yaw_IMUOF]  =QuatToEuler(POSE_IMUOF(:,4:7));
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
title('Orientation of Optical Flow')
plot(OpticalFlowTime(OpticalFlowstart:end),unwrap(roll_OF(OpticalFlowstart:end)));
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),unwrap(roll_IMUOF(OpticalFlowstart:end)),'r');
plot(IMUTime,roll,'k')
xlabel('seconds'); ylabel('roll');
legend('5DoF','IMU and 3DoF','Ground truth');
grid on

subplot(3,1,2);

plot(OpticalFlowTime(OpticalFlowstart:end),pitch_OF(OpticalFlowstart:end));
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),pitch_IMUOF(OpticalFlowstart:end),'r');
plot(IMUTime,pitch,'k')
xlabel('seconds'); ylabel('pitch');
legend('5DoF','IMU and 3DoF','Ground truth');
grid on

subplot(3,1,3);

plot(OpticalFlowTime(OpticalFlowstart:end),yaw_OF(OpticalFlowstart:end));
hold on;
plot(OpticalFlowTime(OpticalFlowstart:end),yaw_IMUOF(OpticalFlowstart:end),'r');
plot(IMUTime,unwrap(yaw),'k')
xlabel('seconds'); ylabel('yaw');
legend('5DoF','IMU and 3DoF','Ground truth');
grid on















