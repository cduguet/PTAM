%REad Vicon File for offline processing

function ReadPoseVicon(flowfilename)

% flowfilename = '20120913_1545_PoseVicon.txt';

% open the file
file_handle = fopen(flowfilename,'r');

%read the header and move the pointer(parameter Headerfile of textscan not working OK)
header = textscan(file_handle,'%s',1,'EndofLine','\n');

[data_read]=textscan(file_handle, '%d %d %d %d %d %f%f%f%f%f%f%f', 'Delimiter',',','EndofLine','\n');

ViconSysTime = [data_read{2:3}];
ViconTime = [data_read{4:5}];
ViconPose=[data_read{6:end}];

save(sprintf('%s.mat',flowfilename),'ViconSysTime','ViconTime','ViconPose');

end