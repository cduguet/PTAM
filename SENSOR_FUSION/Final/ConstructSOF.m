% Construct output of the Optical Flow measurements: 
% How the output should look like 

function [SOF,SOFrotvel,SOFTime,SOFSysTime] = ConstructSOF(ViconTime, ViconSysTime, ViconPose, OFTime, OFSysTime, q_vic_c, p_vic_c)

%assuming small angle assumption interpolate linearly (also quaternions)
PoseInterp = interp1(ViconTime, ViconPose, OFTime);

% calculate velocity w.r.t world frame

vel = diff(PoseInterp(:,1:3)) ./ repmat(diff(OFTime),1,3);

timedifference = diff(OFTime); % sec
timedifference = repmat(timedifference,1,3);

%Convert relative orientation to rot mat
C_q_vic_c = QuatToRotMat(q_vic_c);


%%Calculate Rotation Velocity
% Take the differential quaternion
invPose = invquat(PoseInterp(1:end,4:7)); %inverse quaternion
SOFrotvel = quatmult(PoseInterp(2:end,4:7),invPose(1:end-1,:));

%Take the angular velocity from the differential quaternion.
SOFrotvel = 2*SOFrotvel(1:end,1:3)./timedifference(1:end,:);

%rotate the velocity from the vicon local frame to the IMU frame
SOFrotvel = SOFrotvel * C_q_vic_c';


%Everything: gravity, differentiation, rotation, distance effect
% Convert this world velocity into camera velocity w.r.t. its local frame
for i=1:size(vel,1)
C_q_w_vic = QuatToRotMat(PoseInterp(i+1,4:7));        
vel(i,:) =  (C_q_vic_c * ( C_q_w_vic * vel(i,:)') + C_q_vic_c * (cross(SOFrotvel(i,:)',p_vic_c) ))';
end

SOF=vel;
SOFTime = OFTime(2:end);
% SOFSysTime = interp1(ViconTime,ViconSysTime,OFTime);
SOFSysTime = OFSysTime(2:end);
nan_index = find(isnan(SOF(:,1)));


SOFTime = SOFTime;
SOFSysTime = SOFSysTime;
SOF(nan_index,:)=0;

end 
