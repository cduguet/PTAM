%read and adjust data between vicon and optical flow 

clear all
close all

% function time = readPose(filename)
% read files
% Viconfilename = 'PoseVicon220120315_1754.txt';
% OFfilename = 'log2Dof_20120315175406.txt';
Viconfilename = 'PoseVicon220120315_1756.txt';
OFfilename = 'log2Dof_20120315175656.txt';


Viconfile=fopen(Viconfilename);
OFfile= fopen(OFfilename);

%scan the data 
Vicon=textscan(Viconfile,repmat('%f',1,7),'headerlines',1,'delimiter',',');
OF=textscan(OFfile,repmat('%f',1,3),'headerlines',1,'delimiter',',');




%%%FIXED MANUALLY TO ADJUST DATA TIME AND END (EXPERIMENT WAS NOT
%%%SYNCHRONIZED)
OF{3} = -OF{3};


figure(1)
subplot(2,2,1); plot(Vicon{2}(1:floor(end*0.9)),Vicon{3}(1:floor(end*0.9))); title('Vicon Path')
subplot(2,2,2); plot(OF{2}, OF{3}); title('Optical Flow Path');

subplot(2,2,3); plot(Vicon{2});
subplot(2,2,4); plot(OF{2});

figure(2)
plot(Vicon{4});
% curve corners for OF
%for data 1 
% Viconcorners = [ 686 1876 3216 3771 ];
% OFcorners = [ 32 120 220 261];
% for sdata 2
% Viconcorners = [1316 2298 3346 4060  ];
% OFcorners = [35  94 164 227 ]; %DOESNT WORK
% for sdata 3
Viconcorners = [2039 3210 4107 4911];
OFcorners = [26 113 187 246]; %DOESNT WORK

1

% 775 6395 

diff(Viconcorners)./ diff(OFcorners);

Mbegin = floor(Viconcorners(1) - 13.5 *(OFcorners(1)-1))
Mend = ceil(Viconcorners(end) + 13.5 * (length(OF{1}) -OFcorners(end)))

Mbegin = 775; Mend = 6395;
M = cell(size(Vicon));
for i = 1:length(Vicon)
    M{i} = M{i}(Mbegin:Mend);
end

%compare the form of the curves 
figure(1)
subplot(2,3,1); plot(Vicon{2},Vicon{3}); title('Vicon Path')
subplot(2,3,2); plot(M{2},M{3}); title('Vicon Path')
subplot(2,3,3); plot(OF{2}, OF{3}); title('Optical Flow Path');
subplot(2,3,4); plot(Vicon{4});
subplot(2,3,5); plot(M{4})
subplot(2,3,6); plot(OF{2});








% Interpolate Path of Vicon
 OF{1}=linspace(M{1}(1),M{1}(end),length(OF{1}))';
% 
for i=2:length(M)
    M{i}= spline(M{1},M{i},OF{1});
end
M{1} =OF{1};
M{2} = M{2} - M{2}(1);
M{3} = M{3} - M{3}(1);
save



% once we have the data, we make a particle filter. 

