% Continuous 8 point algorithm, for 3 Degrees of Freedom:
% This function implements the 8 point algorithm, but with the aditional
% constraint that there is no rotation velocity ,i.e. the orientations
% stays constant. To use this function you have to make sure that the
% Optical Flow is unrotated with the use of the IMU Data


function [vel0]= Continuous8point2DoF(x,y,z,u,v,w)

    %Step 1: Define matrix A constructed with the Optical Flow data
  
A = [w.*y-v.*z     u.*z-w.*x    v.*x-u.*y ];
 
    %Determine the SVD of A, and solve Ae=0
    [U,~,V]=svd(A);
    
    % positive  depth constraint
    if [mean(u) mean(v)]*V(1:2,end) >0
        V=-V; U=-U;
    end
    e = V(:,end);
    e = e/norm(e(1:3));

    %recover the velocity vector and symmetric matrix  s 
    vel0 = e(1:3);
    