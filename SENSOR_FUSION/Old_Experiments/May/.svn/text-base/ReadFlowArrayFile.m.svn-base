%Motion calculation in MAtlab
%Author:    Cristian Duguet
% Reads the files saved in ROS by the function flowxy of the package
%This functions is intended to make a .mat file out of the logfile, in
%order to try the optical flow methods in Matlab in an easier way, before
%implementing in C++ 

%temporal: using a filename  instead of a function argument , i.e. we are
%making first a script, which will be transformed into a matlab function
%afterwards


function  ReadFlowArrayFile(flowfilename)

% flowfilename = '20120509_1920_Flow.txt';


% open the file
file_handle = fopen(flowfilename);

% scan the header in lines
[header_text] = textscan(file_handle, '%s',4, 'Delimiter', '\n');
%read FPS
[data_read] = textscan(header_text{1}{2}, '%s %d', 'delimiter',',');
FPS = data_read{2};
%read width and height of frame
[data_read] = textscan(header_text{1}{3}, '%s %d %s %d', 'delimiter',',');
width = data_read{2};
height = data_read{4};
%read maximum number of features detectable
[data_read] = textscan(header_text{1}{4}, '%s %d', 'delimiter',',');
max_feat = data_read{2};

header =textscan(file_handle,'%s',1);

%REad actual data
[data_read]=textscan(file_handle, cat(2,'%d %d ', repmat('%d',1,4*max_feat)),'Delimiter',',');

OpticalFlowTime = [data_read{1}];
matched_feat = [data_read{2}];
OpticalFlow=[data_read{3:end}];


save(sprintf('%s.mat',flowfilename),'OpticalFlowTime', 'max_feat','OpticalFlow','FPS','width','height','matched_feat');

end