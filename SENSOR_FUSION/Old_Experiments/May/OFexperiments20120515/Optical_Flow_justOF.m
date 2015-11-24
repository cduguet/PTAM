%Human and Pedestrian Nagivation homework
%Localization taskusing particle filters 
%Cristian Duguet 

clc
close all
clear all

%Gravitational Constant
G =  9.80665;


experiment_date = '20120523_1724';


%% Parameters:
%in mm
sigmaIMUacc = [0.01 0.01 0.01]; % in (mm/s2)^2

sigmaIMUrotvel = deg2rad([2 2 2]); % in rad^2

%Focal Length (mm):
f = 6; %mm 


% IMU Update Rate:
IMU_rate = 50; %Hz

%Pixel size 
pixelsize=3.6/960;

%% Read files


ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));

load(strcat(experiment_date,'_Flow.txt.mat'));
load(strcat(experiment_date,'_PoseVicon.txt.mat'));


%We noticed that we have asynchronized video and IMU stream. So what we do
%is to make an algorithm which gets this data asynchronized???
    

%% IMU Construction


% --- Acceleration: 
% Build an IMU sensor based on the Vicon position
% FIXME: IMU data is w.r. to the world's fixed frame. Thus it is not totally
% fidel

timedifference = repmat(double(diff(ViconTime)),1,3) /1e9; % sec
% Calculate velocity 
vel = diff(ViconPose(:,1:3)) ./timedifference; %mm/sec
IMUacc = diff(vel) ./timedifference(2:end,:) ; %mm/sec2

% Add Gravity
IMUacc = IMUacc + repmat([0 0 G*1000],size(IMUacc,1),1);


% --- Rotacion velocity: 
% The IMU is delivering the results in quaternions

% Data Correction: We need to rotate the ViconPose by [1 0 0 0] (rot 90deg in x-axis)
ViconPose(4:7) = quatmult([1 0 0 0],ViconPose(4:7));


%FIXME: rotation velocity is given just as difference of quaternions, but
%not as rotation velocity
% Take the differential quaternion 
invPose = invquat(ViconPose(1:end-1,4:7));
IMUrotvel = quatmult(ViconPose(2:end,4:7),invPose);

IMUTime_instants = size(IMUacc,1);


% ---- Noise: 

% Translational Noise 
IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigmaIMUacc,IMUTime_instants,1));

% Rotational Noise

% Create Noise Quaternion
% Using the small angle assumtion, we can model the noise quaternion as a 
% q_noise = (noise_x noise_y noise_z 1); assuming noise_i ~0

NoiseQuat = [ random('Normal', zeros(IMUTime_instants,3)    , repmat([sigmaIMUrotvel],IMUTime_instants,1))        ones(IMUTime_instants,1)];

IMUrotvel = quatmult(IMUrotvel(2:end,:),NoiseQuat);


%% Time adjustments 

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

IMUTime = ViconTime(3:end);

% the time in which both are transmiting
time_zero =  min(OpticalFlowTime(1), IMUTime(1));

OpticalFlowTime = double(OpticalFlowTime-time_zero) * 1e-9;
IMUTime = double(IMUTime - time_zero) * 1e-9;


% ---------- Test ---------------------------------------------------------
% figure(1)
% plot(ViconPose(:,1),ViconPose(:,2))
% 
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


%% Odometry 

% FIXME: I don't know yet but the sample rate of the camera does not seem
% to be 15FPS 

% display(sprintf('Mean rate of flor information: %g arrays per second \n',1e9/(mean(diff(OpticalFlowTime)))));


% Here ends the data preparation and begins the Odometry algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Odometry and Predicion 

% * Considerations:

% Given that the flow velocities (u,v) are not expressed related to time, but
% rather as differences, the angular velocities and the translational
% velocities must be as well
 
% let's make the time steps equivalent to the biggest one, which is in this
% case the Camera frame rate. The IMU data will be integrated until its
% sample right after the new OpticalFlow sample timestamp. 

% Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
X = [0 0 1000 0 0 0 1]; 

%% First iterations
% We are not going to use the Optical Flow data which comes after we
% started receiving the IMU data. 

OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);

% -------------------------- Optical Flow----------------------------------
% --------------------------    IMU     -----------------------------------
 
%we don't use IMUacc(0,:)

% Initial velocity: Very important!!! Otherwise the position error will be
% propagating
% ODO_IMUvel = vel(2,:);  %[0 0 0];
% ODO_IMUvel = IMUacc(1,:) .* double(ViconTime(3) - ViconTime(2)) * 1e-9;
ODO_IMUvel = [0 0 0];


% start at odometry movement (0,0,0) and rotation (0,0,0,1) (neutral
% quaternion)
% ODO_IMUpos = ODO_IMUvel.* double(ViconTime(3) - ViconTime(2)) * 1e-9
ODO_IMUpos = [0 0 0];
ODO_IMUrot = [ 0 0 0 1]; % no rotation

% --- Test
% fprintf('Debugging: \t Step \t Time Difference  \t Acc \t Vel \n' );

for j = 2: IMUsecondit
    % Translation
    delta_t_IMU = IMUTime(j) - IMUTime(j-1);
    
    ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU; 
    ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
    
    % Rotation
    ODO_IMUrot = quatmult(IMUrotvel(j,:),ODO_IMUrot);
end

% -------------------------------------------------------------------------
% ------------ Test
testpos_OF3D = [0 0];
testpos_OF2D = [0 0];
testpos_IMU = [0 0];
% 
figure(6)
% 
actualpos3  = testpos_IMU;
testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
% 
plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
hold on
% -------------------------------------------------------------------------


% %  UPDATE POSE!!!!! 
X = [ X(1:3) + ODO_IMUpos quatmult(ODO_IMUrot,X(4:7))];

%% Second iteration and So on
index_end = IMUsecondit;
for i = OpticalFlowstart+1: length(OpticalFlowTime)
    
    %% Calculate the IMU movement : 
    
    index_start = index_end; %use the end index of last iteration
    index_end = find(IMUTime > OpticalFlowTime(i),1); % FIXME: can crash if theres no IMU sample after the last OF sample 
    
    ODO_IMUpos = [0 0 0];
    ODO_IMUrot = [ 0 0 0 1]; % no rotation 
    
    for j = index_start+1: index_end    
        % Translation
        delta_t_IMU = IMUTime(j) - IMUTime(j-1);
        ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU;
        ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
        
        % Rotation
        ODO_IMUrot = quatmult(IMUrotvel(j,:),ODO_IMUrot);
        %        --- Test
        %     if j<23
        %     fprintf('\t %u | %d | %d  %d | %d %d  | %d %d | \n\n ',j,delta_t_IMU,IMUacc(j,1:2),ODO_IMUvel(1:2),ODO_IMUpos(1:2)); end
    end
    
    
    
    %% Optical Flow
    
    % Set the pixel position referenced to the center
    x= double(OpticalFlow_x(i,1:matched_feat(i)))'-Image_width/2-0.5;
    y= double(OpticalFlow_y(i,1:matched_feat(i)))'-Image_height/2-0.5;
    
    
    % Construct the Z vector
    % Determine the distance of every pixel / feature and save it to Z
    Z = zeros(matched_feat(i),1);
    
    for j = 1: matched_feat(i)
        
        % we first construct a rotation quaternion representing the
        % coordinate conversion from a picel to the central camera frame
        % using the pixel position and the focal length
        
        r=sqrt(x(j).^2+y(j).^2);
        angle_camera = tan(pixelsize*r/f);
        
        q_frompixeltocamera= [sin(angle_camera)/r*[y(j) -x(j) 0] cos(angle_camera)];
        
        % Then we construct the projection quaternion from the pixel to
        % world. Here X(4:7) is considered as the quaternion from camera to
        % world
        q_proj = quatmult(X(4:7),q_frompixeltocamera);
        
        % A vector pointing in the -z direction in the pixel frame
        % represented in the worlds frame is 
        v_projected = quatmult(q_proj,[0 0 -1 0]);
        v_projected = quatmult(v_projected,invquat(q_proj));

        
        % Then, the angle difference to the world's vertical is 
        % (we consider that the vectors and quaternions are normalized)
        angle_world = acos([0 0 -1]*v_projected(1:3)');
        
        % Then, the Z component of the pixel j is 
        Z(j)  = X(3)/cos(angle_world) * cos(angle_camera);
    end

    
    
    %MAKE Z!!!! 
    
    %      [u;v] = [-f/Z 0  x/Z  x*y/f   -f-x^2/f y;
    %                0 -f/Z y/Z f+y^2/f  -x*y/f   x] * posvector;
    
    U = double([OpticalFlow_u(i,1:matched_feat(i)) OpticalFlow_v(i,1:matched_feat(i))]');                    

    P = [[  -f./Z                zeros(matched_feat(i),1)    x./Z      x.*y./f         -f-x.^2./f      y];
        [ zeros(matched_feat(i),1)        -f./Z             y./Z     f+y.^2./f          -x.*y./f      -x]];
        
    %Least squares
    pos = P\U;
    
    
    %Convert the angular velocity 
    ODO_OpticalFlowpos = pos(1:3)';
    ODO_OpticalFlowrot = [pos(4:6)'/2 1]; % FIXME!!!!!! USING SMALL ANGLE ASSUMPTION      
   % The Optical Flow image motion equations give us the angular
   % velocities, but these ara actually not measured in rad/seg, because
   % the optical flow velocity isnt measured in pix/seg either!!! What it
   % gives us isjus ta difference. In which order do these rotations
   % occur??? It looks like the Image velocity equations have already taken
   % the small angle asumption, and that is why it looks like the order
   % doesnt matter . Still... not sure about it. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%%  ------------ Test -
        
%     OpticalFlow2D

%     Assume That the camera movement is uniform in the camera. Orientation is fixed. height is fixed. 

    ODO_OpticalFlow2Dpos = [mean(double(OpticalFlow_u(i,1:matched_feat(i))))*X(3)/f  mean(double(OpticalFlow_u(i,1:matched_feat(i))))*X(3)/f];
    
    %tests movements in 2D
%     figure(1)
%     actualpos = testpos_OF3D;
%     testpos_OF3D =  testpos_OF3D + ODO_OpticalFlowpos(1:2);
%     
%     plot([actualpos(1) testpos_OF3D(1) ],[actualpos(2) testpos_OF3D(2)])
%     hold on 
    
% test remove outliers 


    figure(4)
    
    actualpos2 = testpos_OF2D;
    testpos_OF2D =  testpos_OF2D + ODO_OpticalFlow2Dpos(1:2);
    plot([actualpos2(1) testpos_OF2D(1) ],[actualpos2(2) testpos_OF2D(2)])
    hold on
    
    figure(5)
    plot([ x' ;  x'+double(OpticalFlow_u(i,1:matched_feat(i)))],[y' ;  y'+double(OpticalFlow_v(i,1:matched_feat(i)))])

%     figure(5)
%     plot(ViconPose(:,1),ViconPose(:,2));

figure(6)

actualpos3  = testpos_IMU;
testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
% if abs(testpos_IMU - k(j,1:2)) > 1e-5
%     error('the do not coincide');
% else
%     fprintf('the IMU integrations coincide!! :):) \n');
% end
plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
hold on
% ---------------------------------------------------------------------
%% Kalman Filter
end
