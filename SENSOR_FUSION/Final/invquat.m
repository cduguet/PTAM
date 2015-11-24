% Inverse of a Quaterion 
%   by Cristian Duguet
% Info: The input Q can be a quaternion in the form [x y z w] , i.e a row
% vector or an array of N quaternions in a Matrix of size Nx4
function Qinv =invquat(Q)
    N = size(Q,1);
    Qinv = zeros(N,4);
    
    Qinv = [-Q(:,1:3) Q(:,4)];
    
return