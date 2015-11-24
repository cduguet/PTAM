%Create and try quaternion stuff

clc
clear all 


%  camera rotation 
vectorrot = [3 -2 4]

angle = deg2rad(30)

% q_c^w : from camera coord frame to world coord frame
q1 = [vectorrot/norm(vectorrot)*sin(angle/2) cos(angle/2)]

% pixel projection angle: rotation angle relative to the camera centered
% coordination system 

vectorrot = [1 1 0]
angle = deg2rad(0)

% q_p^c : from pixel projection to camera 
q2 = [vectorrot/norm(vectorrot)*sin(angle/2) cos(angle/2)]

%then we need q_p^w 


% now we multiplicate these two rotations
q3 = quatmult(q1,q2);


% now we try to find out what is the angle w.r.t. the vertical 

angle1 = 2*atan2(norm(q3(1:3)),q3(4)); % i suspect this angle is the angle considering the z rotation it is bigger than what i am looking for 

angle1 = rad2deg(angle1)

%with this we want to know the angle between a vector (0 0 -1) in the pixel
%cf and a vector (0 0 -1) in the world c.f.


verticalproj = quatmult(q3,[0 0 -1 0]);
verticalproj = quatmult(verticalproj,invquat(q3));


%then the angle difference is 

angle2 = acos([0 0 -1]*verticalproj(1:3)');
angle2 = rad2deg(angle2)
    



