%ReadIMU 

function ReadIMU(imufilename)

% imufilename = '20120913_1545_IMU.txt';

% open the file
file_handle = fopen(imufilename,'r');


% scan the header in lines
[header_text] = textscan(file_handle, '%s',4, 'Delimiter', '\n','EndofLine','\n');
%read Angular Velocity Covariance
[data_read] = textscan(header_text{1}{2}, '%s %f%f%f%f%f%f%f%f%f', 'delimiter',',');
AngularVelocityCovariance = reshape([data_read{2:end}],3,3)';
%read Acceleration Covariance
[data_read] = textscan(header_text{1}{3}, '%s %f%f%f%f%f%f%f%f%f', 'delimiter',',');
AccelerationCovariance = reshape([data_read{2:end}],3,3)';

%read the header and move the pointer(parameter Headerfile of textscan not working OK)
[data_read] = textscan(file_handle,'%d %d %d %d %d %f%f%f%f%f%f%f%f%f%f','Delimiter',',','EndofLine','\n');

IMUSysTime = [data_read{2:3}];
IMUTime = [data_read{4:5}];
IMUrot = [data_read{6:9}];
IMUrotvel = [data_read{10:12}];
IMUacc = [data_read{13:15}];

save(sprintf('%s.mat',imufilename),'IMUSysTime','IMUTime','IMUrot','IMUrotvel','IMUacc','AngularVelocityCovariance', 'AccelerationCovariance');

end