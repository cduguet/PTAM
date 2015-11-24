%Motion calculation in MAtlab
%Author:    Cristian Duguet
% Reads the files saved in ROS by the function flowxy of the package
% This functions saves the data saved by the OpticalFlow logger program,
% particulartly the camera information. The camera_info file is important
% because it includes the time stamps of the camera, with some other
% additional information 


function  ReadCameraInfo(infofilename)
% flowfilename = '20120913_1545_Flow.txt';

% open the file
file_handle = fopen(infofilename);

%Read actual data
[data_read]=textscan(file_handle, '%u %d %d %d %d','Delimiter',',','HeaderLines',1,'EndofLine','\n');

CameraInfoSysTime = [data_read{2:3}];
CameraInfoTime = [data_read{4:5}];

save(sprintf('%s.mat',infofilename),'CameraInfoSysTime','CameraInfoTime');
end