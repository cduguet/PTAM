%Rotation order roll pitch yaw

function [roll,pitch yaw] = QuatToEuler(q)
roll  =zeros(size(q,1),1);
pitch  =zeros(size(q,1),1);
yaw =zeros(size(q,1),1);

for i=1:size(q,1)
    
roll(i) = atan2( 2 * (q(i,4)*q(i,1) + q(i,2)*q(i,3)), 1- 2 *(q(i,4)^2 + q(i,2)^2));
pitch(i) = asin(2*(q(i,4)*q(i,2) - q(i,3)*q(i,1)));
yaw(i) = atan2( 2 * (q(i,4)*q(i,3) + q(i,2)*q(i,1)), 1- 2 *(q(i,3)^2 + q(i,2)^2));
end
