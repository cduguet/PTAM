 
function [SIMUacc, SIMUrotvel,SIMUTime, SIMUSysTime] = Construct_SIMU(ViconTime,ViconSysTime,ViconPose,q_vicon_i,p_vicon_i)
% Gravitational Constant
G =  9.80665;

% Build an IMU sensor based on the Vicon Pose
timedifference = diff(ViconTime); % sec
timedifference = repmat(timedifference,1,3);



C_q_vic_i = QuatToRotMat(q_vicon_i);


% --- Rotation velocity:
% The IMU is delivering the results in angular velocity

% Data Correction: We need to rotate the ViconPose by
% FIXME: Just for this experimen
%Rotation -90 in Z axis then 90 in X axis: [.5 -.5 -.5 .5]
% ViconPose(:,4:7) = quatmult(repmat([.5 -.5 -.5 .5],size(ViconPose,1),1),ViconPose(:,4:7));
%FIXME: maybe its better to set the camera orientation as looking downwards
%([-1 0 0 0])

% Take the differential quaternion
invPose = invquat(ViconPose(1:end,4:7)); %inverse quaternion
SIMUrotvel = quatmult(ViconPose(2:end,4:7),invPose(1:end-1,:));
% IMUrotvel = quat_normalize(IMUrotvel);

%FIXME: Vicon sometimes deliver bad irregular rotations (due to jumps). Supress these
%making it equal to their previous rotation velocity
for i=1:size(SIMUrotvel,1)
    if abs(SIMUrotvel(i,4) - 1) > 0.1
        SIMUrotvel(i,:) = [0 0 0 1];
        ViconPose(i+1,4:7) = ViconPose(i,4:7);
        invPose(i,:) = invPose(i-1,:);
        error('error: Vicon has jumps \n');
    end
end

%Take the angular velocity from the differential quaternion.
SIMUrotvel = 2*SIMUrotvel(1:end,1:3)./timedifference(1:end,:);

%rotate the velocity from the vicon local frame to the IMU frame
SIMUrotvel = SIMUrotvel * C_q_vic_i';

% --- Acceleration -----:

% Calculate velocity
vel = diff(ViconPose(:,1:3))./timedifference; %m/sec


SIMUacc = zeros(size(vel) - [1 0]);

a_g = [0 0 -G]';

%Everything: gravity, differentiation, rotation, distance effect
% Convert this world velocity into the imu velocity w.r.t. its local frame
for i=1:size(vel,1)-1
C_q_w_vic = QuatToRotMat(ViconPose(i,4:7));        
delta_t  = timedifference(i+1,1);    
% SIMUacc(i,:) =   ( C_q_w_vic * ((vel(i+1,:)-vel(i,:))' + a_g));
SIMUacc(i,:) =  C_q_vic_i * ( C_q_w_vic * ((vel(i+1,:)-vel(i,:))'/delta_t  + a_g));
% SIMUacc(i,:) = SIMUacc(i,:) - (C_q_vic_i * (cross(SIMUrotvel(i,:)',C_q_w_vic*vel(i+1,:)'))/delta_t)';
% SIMUacc(i,:) = SIMUacc(i,:) + (C_q_vic_i * (cross(SIMUrotvel(i+1,:)'-SIMUrotvel(i,:)',p_vicon_i) ) /delta_t )';
end
SIMUrotvel = SIMUrotvel(2:end,:);
% 
% figure;
% plot(vel(:,1:2))
% 
% figure;plot(cumsum(SIMUacc(:,1:2).*timedifference(2:end,1:2)))
% 
% figure;plot(cumsum(diff(vel(:,1:2))))

% %%-- debug
% SIMUacc_no_g = SIMUacc;
% SIMUacc_no_g = quatmult([SIMUacc_no_g zeros(size(SIMUacc_no_g,1),1)],invPose(3:end,:));
% SIMUacc_no_g = quatmult(ViconPose(3:end,4:7),SIMUacc_no_g);
% SIMUacc_no_g = SIMUacc_no_g(:,1:3);
% %%- debug

% ---- Noise:
%FIXME; No bias considered
% Translational Noise
% SIMUacc = SIMUacc + random('Normal',zeros(SIMUTime_instants,3), repmat(diag(sigmaIMUacc)',SIMUTime_instants,1));

% Rotational Noise
% SIMUrotvel = SIMUrotvel + random('Normal', zeros(SIMUTime_instants,3), repmat(diag(sigmaIMUrotvel)',SIMUTime_instants,1));

SIMUTime  = ViconTime(3:end);
SIMUSysTime = ViconSysTime(3:end);
end