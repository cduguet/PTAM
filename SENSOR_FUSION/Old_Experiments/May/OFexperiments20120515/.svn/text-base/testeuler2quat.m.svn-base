clc, close all, clear all
format 

dx=160; dy=100; dz=120;               % euler angles in degree  

D2R=pi/180;                           % for converting degree to radians
rx=dx.*D2R; ry=dy.*D2R;   rz=dz.*D2R; % euler angles in radinas  

% % x i.e. roll
% % y i.e. pitch
% % z i.e. yaw
% % Quaterion is expressed as Q=[w x y z]
% % w is scalar component of quaternion
% % x,y,z are vector components of quaternion. For example q=w+xi+yj+zk

% % euler2quat converts euler angles to quaternion for six basic sequence
% % of rotations.
% % Allowed rotation sequences: xyz, xzy, yxz, yzx, zxy, zyx 

[w,x,y,z]=euler2quat(rx,ry,rz,'xyz');  %order of rotation is 'xyz'
q=[w x y z]

[w,x,y,z]=euler2quat(rx,ry,rz,'xzy');  %order of rotation is 'xzy'
q=[w x y z]

[w,x,y,z]=euler2quat(rx,ry,rz,'yxz');  %order of rotation is 'yxz'
q=[w x y z]

[w,x,y,z]=euler2quat(rx,ry,rz,'yzx');  %order of rotation is 'yzx'
q=[w x y z]

[w,x,y,z]=euler2quat(rx,ry,rz,'zxy');  %order of rotation is 'zxy'
q=[w x y z]

[w,x,y,z]=euler2quat(rx,ry,rz,'zyx'); %order of rotation is 'zyx'
q=[w x y z]



% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------

% Technical Reference: Ken Shoemake, "Animating Rotations with Quaternion Curves"


