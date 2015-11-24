% Conver from quaternion to Rotation MAtrix 
% form of quaternion must be [x y z w]' (horizontal or vertical vector)

function C = QuatToRotMat(q)

if length(q) ~=4
    error('Quaternion has not the form [x y z w]');
    return
end

crossprodmat = VectorToCrossMat(q(1:3));
            
C = eye(3) - 2*q(4)*crossprodmat + 2*crossprodmat^2;
end