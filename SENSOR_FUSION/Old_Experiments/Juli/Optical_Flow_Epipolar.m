%Human and Pedestrian Nagivation homework
%Localization taskusing particle filters 
%Cristian Duguet 

clc
close all
clear all

%Gravitational Constant
G =  9.80665;


experiment_date = '20120704_1913';


%% Parameters:
%in mm
sigmaIMUacc = [0.01 0.01 0.01]; % in (mm/s2)^2
sigmaIMUrotvel = deg2rad([2 2 2]); % in rad^2

%Focal Length (mm):
f = 6; %mm 

% IMU Update Rate:
IMU_rate = 50; %Hz

%Scaling constant to be changed 
scaling = 1/40; 

MAX_NOFLOW = 100;

%% Read files

%FIXME: PoseVicon logger is not saving the w componnent!!!! 
try
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
catch ME
    ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
    ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));
    
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
end

delta = diff(ViconTime,1);
ViconFPS = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

delta = diff(OpticalFlowTime);
OpticalFlowFPS = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));


%We noticed that we have asynchronized video and IMU stream. So what we do
%is to make an algorithm which gets this data asynchronized???
 
%FIXME: number of columns are less.... 
%max_features is one less i think 
%print 0.0 instead of 0 
% OpticalFlow_u = OpticalFlow_u(:,1:end-2);
% OpticalFlow_v = OpticalFlow_v(:,1:end-2);
% OpticalFlow_x = OpticalFlow_x(:,1:end-2);
% OpticalFlow_y = OpticalFlow_y(:,1:end-2);

max_feat=398;

matched_feat(find(matched_feat>max_feat)) = max_feat;

%Coverting ground truth into mm
ViconPose(:,1:3) = ViconPose(:,1:3)*1000;

%Pixel size 
pixel_size=3.6/Image_width;

noflow_counter =0;

%% IMU Construction


% --- Acceleration: 
% Build an IMU sensor based on the Vicon position
% FIXME: IMU data is w.r. to the world's fixed frame. Thus it is not totally
% fidel

timedifference = double(diff(ViconTime(1,:))) + double(diff(ViconTime(:,2)))/1e9; % sec
timedifference = repmat(timedifference,1,3);
% Calculate velocity 

vel = diff(ViconPose(:,1:3)) ./timedifference; %mm/sec
IMUacc = diff(vel) ./timedifference(2:end,:) ; %mm/sec2

% Add Gravity
IMUacc = IMUacc + repmat([0 0 G*1000],size(IMUacc,1),1);


% --- Rotacion velocity: 
% The IMU is delivering the results in quaternions

% Data Correction: We need to rotate the ViconPose by [1 0 0 0] (rot 90deg in x-axis)
ViconPose(4:7) = quatmult([1 0 0 0],ViconPose(4:7));


%FIXME: rotation velocity is given just as difference of quaternions, but
%not as rotation velocity
% Take the differential quaternion 
invPose = invquat(ViconPose(1:end-1,4:7));
IMUrotvel = quatmult(ViconPose(2:end,4:7),invPose);

IMUTime_instants = size(IMUacc,1);


% ---- Noise: 

% Translational Noise 
IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigmaIMUacc,IMUTime_instants,1));

% Rotational Noise

% Create Noise Quaternion
% Using the small angle assumtion, we can model the noise quaternion as a 
% q_noise = (noise_x noise_y noise_z 1); assuming noise_i ~0

NoiseQuat = [ random('Normal', zeros(IMUTime_instants,3)    , repmat([sigmaIMUrotvel],IMUTime_instants,1))        ones(IMUTime_instants,1)];

IMUrotvel = quatmult(IMUrotvel(2:end,:),NoiseQuat);


%% Time adjustments 

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

IMUTime = ViconTime(3:end);

% the time in which both are transmiting
time_zero =  min(OpticalFlowTime(1), IMUTime(1));

OpticalFlowTime = double(OpticalFlowTime-time_zero) * 1e-9;
IMUTime = double(IMUTime - time_zero) * 1e-9;


% ---------- Test ---------------------------------------------------------
% figure(1)
% plot(ViconPose(:,1),ViconPose(:,2))

% Test: to test how it would be different if we know the initial
% velocity
% h = cumsum((cumsum(IMUacc .*timedifference(2:end,:),1)+repmat(vel(1,:),size(IMUacc,1),1))  .*timedifference(2:end,:), 1);
% figure(2)
% plot(h(:,1),h(:,2))
% 
% kk = cumsum(IMUacc.*timedifference(2:end,:) ,1);
% k = cumsum( kk.*timedifference(2:end,:) ,1) ;
% figure(3)
% plot(k(:,1),k(:,2))
% 
% 
% fprintf('Debugging: \t Step \t Time Difference  \t Acc \t Vel \n' );
% for i=2:23
% fprintf('\t %u | %d | %d %d | %d %d | %d %d | \n\n ',i,timedifference(i+1,1),IMUacc(i,1:2), kk(i,1:2), k(i,1:2));
% end


%% Odometry 

% FIXME: I don't know yet but the sample rate of the camera does not seem
% to be 15FPS 

% display(sprintf('Mean rate of flor information: %g arrays per second \n',1e9/(mean(diff(OpticalFlowTime)))));


% Here ends the data preparation and begins the Odometry algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Odometry and Predicion 

% * Considerations:

% Given that the flow velocities (u,v) are not expressed related to time, but
% rather as differences, the angular velocities and the translational
% velocities must be as well
 
% let's make the time steps equivalent to the biggest one, which is in this
% case the Camera frame rate. The IMU data will be integrated until its
% sample right after the new OpticalFlow sample timestamp. 

% Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
%FIXME: We are taking the ground truth initial position, for comparison
X = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1); 

%% First iterations
% We are not going to use the Optical Flow data which comes after we
% started receiving the IMU data. 

OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);

% -------------------------- Optical Flow----------------------------------
% --------------------------    IMU     -----------------------------------
 
%we don't use IMUacc(0,:)

% Initial velocity: Very important!!! Otherwise the position error will be
% propagating
% ODO_IMUvel = vel(2,:);  %[0 0 0];
% ODO_IMUvel = IMUacc(1,:) .* double(ViconTime(3) - ViconTime(2)) * 1e-9;
ODO_IMUvel = [0 0 0];


% start at odometry movement (0,0,0) and rotation (0,0,0,1) (neutral
% quaternion)
% ODO_IMUpos = ODO_IMUvel.* double(ViconTime(3) - ViconTime(2)) * 1e-9
ODO_IMUpos = [0 0 0];
ODO_IMUrot = [ 0 0 0 1]; % no rotation

% --- Test
% fprintf('Debugging: \t Step \t Time Difference  \t Acc \t Vel \n' );

for j = 2: IMUsecondit
    % Translation
    delta_t_IMU = IMUTime(j) - IMUTime(j-1);
    
    ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU; 
    ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
    
    % Rotation
    ODO_IMUrot = quatmult(IMUrotvel(j,:),ODO_IMUrot);
end

% -------------------------------------------------------------------------
% ------------ Test
% testpos_OF3D = [0 0];
% testpos_OF2D = [0 0];
% testpos_IMU = [0 0];
% % 
% figure(6)
% % 
% actualpos3  = testpos_IMU;
% testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
% % 
% plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
% hold on
% -------------------------------------------------------------------------



%% Second iteration and So on : here we include the Optical Flow
max_inliers = zeros(size(matched_feat));


index_end = IMUsecondit;
for i = OpticalFlowstart+1: length(OpticalFlowTime)
    
    %% Calculate the IMU movement :   
    index_start = index_end; %use the end index of last iteration
    index_end = find(IMUTime > OpticalFlowTime(i),1); % FIXME: can crash if theres no IMU sample after the last OF sample 
    
    ODO_IMUpos = [0 0 0];
    ODO_IMUrot = [ 0 0 0 1]; % no rotation 
    
    for j = index_start+1: index_end    
        % Translation
        delta_t_IMU = IMUTime(j) - IMUTime(j-1);
        ODO_IMUvel = ODO_IMUvel + IMUacc(j,:) .* delta_t_IMU;
        ODO_IMUpos = ODO_IMUpos + ODO_IMUvel * delta_t_IMU;
        
        % Rotation
        ODO_IMUrot = quatmult(IMUrotvel(j,:),ODO_IMUrot);
        %        --- Test
        %     if j<23
        %     fprintf('\t %u | %d | %d  %d | %d %d  | %d %d | \n\n ',j,delta_t_IMU,IMUacc(j,1:2),ODO_IMUvel(1:2),ODO_IMUpos(1:2)); end
    end
    
    
    
    %% =====================OPTICAL FLOW===================================
    ama=pixel_size;
    
    % Set the pixel position referenced to the center
    x = (OpticalFlow_x(i,1:matched_feat(i)))'-Image_width/2+0.5;
    y = (OpticalFlow_y(i,1:matched_feat(i)))'-Image_height/2+0.5;
    x = x *(ama/f); y = -y *(ama/f);
    z = ones(size(x));
    
    u = (OpticalFlow_u(i,1:matched_feat(i)))'*(ama/f);
    v = -(OpticalFlow_v(i,1:matched_feat(i)))'*(ama/f);      
    
    %FIXME now they are just being divided by f, but not multiplied by
    %pixelsize
    
    
    %%In case of no flow detected
    if length(u) < 5
        if     noflow_counter > MAX_NOFLOW
            error('Optical Flow information too poor. Need more features.\n'); 
            return
        end
        warning('No flow detected in this image frame.\n');
        noflow_counter  = noflow_counter+1;
        %keep constant movement
        ODO_OpticalFlowpos;
        ODO_OpticalFlowrot;
        X(i,:) = [ X(i-1,1:3)+ODO_OpticalFlowpos quatmult(ODO_OpticalFlowrot,X(i-1,4:7))];
    continue
    end
    
       
    %% ---------------------RANSAC-----------------------------------------
    
    
    
    % First method: Calculating the Homogrpahy Matrix from a random subset
    % of sample_num points of all the matched features in the optical flow.
    % After max_Iter iterations, return the inliers of the most accepted
%     % model (The model which can include the most inliers).
%     max_Iter = 1000;
%     sample_num = 5;
%     RANSAC_Homography_treshold = 1;
%     
%     inliers_idx = RANSAC_Homography(x,y,u,v,max_Iter,sample_num,RANSAC_Homography_treshold);
    
    


    % Second method: Choosing the Consensus of the random samples in a
    % linear sense (without Homography model) 
% 
    max_Iter = 1000;
    sample_num = 50;    
    RANSAC_Linear_treshold = 1;

    inliers_idx = RANSAC_Linear(u,v,max_Iter,sample_num,RANSAC_Linear_treshold);

    
    full_x=x;
    full_y=y;
    full_u=u;
    full_v=v;
    x = x(inliers_idx);
    y = y(inliers_idx);
    z = z(inliers_idx);
    u = u(inliers_idx);
    v = v(inliers_idx);
    
    max_inliers(i) = length(inliers_idx);

    
    %%----------------RANSAC_END-------------------------------------------
    

    %% ---------------EPIPOLAR RECONSTRUCTION------------------------------
    %Step 1: Define matrix A constructed with the Optical Flow data
    %u3 not observable??
    u3=0;
    A = [u3.*y-v.*z     u.*z-u3.*x    v.*x-u.*y   x.^2    2.*x.*y    2.*x.*z   y.^2    2.*y.*z    z.^2];
    
    %case w=(0,0,0)
%     A = [u3.*y-v.*z     u.*z-u3.*x    v.*x-u.*y ];
    %case w=(0,0,0), v3 =0;
%     A = [u3.*y-v.*z     u.*z-u3.*x];
    
    %Determine the SVD of A, and solve Ae=0
    [U,S,V]=svd(A);
    
    % positive  depth constraint
    if [mean(u) mean(v)]*V(1:2,end) >0
        V=-V; U=-U;
    end
    e = V(:,end);
    e = e/norm(e(1:3));

    %recover the velocity vector and symmetric matrix  s 
    vel0 = e(1:3);
    
    s = [[e(4) e(5) e(6)];
         [e(5) e(7) e(8)];
         [e(6) e(8) e(9)]];
%     Eigenvalue decomposition of the symmetric matrix
                
    [V1,lambda] = eig(s);
    lambda = diag(lambda);
    lambda = flipud(lambda);
    V1 = fliplr(V1);
    
    %DEBUG: Dealing with complex eigenvalues 
    % The nearest matrix to a one with complex eingenvalues is the one with
    % real truncated igenvalues [\cite
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/290605]
    if ~isreal(lambda)
     warning('time %d: complex eigenvalues!. Truncating... \n', i);
     return
     lambda = real(lambda);
    end
    
    % eig() function does not give the eigenvalues sorted [\cite 
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/156868]. 
    % We have to do it manually
    if lambda(1) < lambda(2) || lambda(2) < lambda(3)
        [lambda,idx] = sort(lambda,1,'descend');
        V1 = V1(:,idx);
    end
    
    %in case the eigenvalues do not complain with our restrictions:
    lambda(1) = max([lambda(1) 0]);  
    lambda(3) = min([lambda(3) 0]);
    lambda(2) = max([lambda(2) lambda(3)]);
    lambda(2) = min([lambda(2) lambda(1)]);
    
    
    %Projet symmetric matrix onto the symmetric epipolar space S 
    sigma = zeros(3,1);
    sigma(1) = (2*lambda(1) + lambda(2) - lambda(3))/3;
    sigma(2) = (lambda(1) + 2*lambda(2) + lambda(3))/3;
    sigma(3) = (2*lambda(3) + lambda(2) - lambda(1))/3;
    
    %Recover velocity from the symmetric epipolar part
    Lambda = sigma(1) - sigma(3);
    Theta = acos(-sigma(2)/Lambda);
    if Theta > pi || Theta < 0
        warning('Lambda or Theta in symmetric epipolar component are out of bounds');
    end
    
    Ry  = inline(' [[cos(x) 0 sin(x)]; [0  1   0 ]; [-sin(x) 0  cos(x)]] ');
    Rz  = inline(' [[cos(x) -sin(x) 0]; [sin(x) cos(x) 0] ; [0 0 1] ]');
    V = V1 * Ry(Theta/2 - pi/2)';
    U = -V * Ry(Theta);

    SIGMA_Lambda = diag([Lambda Lambda 0]);
    SIGMA_1      = diag([1 1 0]);
    
    %construct matrix with the 4 solutions as columns
    omega = zeros(3,4);
    v_     = zeros(3,4);
    
    % Calculate the 4 different solutions:
    aux = U * Rz(pi/2) * SIGMA_Lambda * U';     omega(:,1) = [aux(8) -aux(7) aux(4)];
    aux = V * Rz(pi/2) * SIGMA_1      * V';     v_(:,1)     = [aux(8) -aux(7) aux(4)];
    
    aux = U * Rz(-pi/2) * SIGMA_Lambda * U';    omega(:,2) = [aux(8) -aux(7) aux(4)];
    aux = V * Rz(-pi/2) * SIGMA_1      * V';    v_(:,2)     = [aux(8) -aux(7) aux(4)];
    
    
    aux = V * Rz(pi/2) * SIGMA_Lambda * V';     omega(:,3) = [aux(8) -aux(7) aux(4)];
    aux = U * Rz(pi/2) * SIGMA_1      * U';     v_(:,3)     = [aux(8) -aux(7) aux(4)];
    
    
    aux = V * Rz(-pi/2) * SIGMA_Lambda * V';    omega(:,4) = [aux(8) -aux(7) aux(4)];
    aux = U * Rz(-pi/2) * SIGMA_1      * U';    v_(:,4)     = [aux(8) -aux(7) aux(4)];
    
    
    % choose the solution with v most similar to vel0
    [aux,idx] = max(vel0'*v_);
    
    omega = omega(:,idx);
    v_ = v_(:,idx);
    %% SCALE RECOVERY
    
    %  crossprod = cross(repmat(omega,1,max_inliers(i)), [x';y';z'],1);
    %  M = [[diag(u)-diag(crossprod(1,:))                  , diag(x), -repmat(vel0(1),max_inliers(i),1)];
    %       [diag(v)-diag(crossprod(2,:))                  , diag(y), -repmat(vel0(2),max_inliers(i),1)];
    %       [diag(u3,max_inliers(i)-1)-diag(crossprod(3,:)), diag(z), -repmat(vel0(3),max_inliers(i),1)]];
    %
    %     [~,~,V2] = svd(M);
    %     lambda_scale = V2(:,end);
    %     scaling = lambda_scale(end);
    %
    
    %% Assign Odometry
    
    % Recover angular motion and velocity 
    %Convert the angular velocity 
    %USE PREVIOUS STATE
    ODO_OpticalFlowpos = vel0' *X(i-1,3) * scaling;
    
    
    
    
    
    ODO_OpticalFlowrot = [omega'/2 1]; % FIXME!!!!!! USING SMALL ANGLE ASSUMPTION     
    
    % %  UPDATE POSE!!!!!   
    if ~isreal(ODO_OpticalFlowrot)
        warning('Oops!... complex numbers in rotation quaternion');
        return
    end
    X(i,:) = [ X(i-1,1:3)+ODO_OpticalFlowpos quatmult(ODO_OpticalFlowrot,X(i-1,4:7))];

    

   % The Optical Flow image motion equations give us the angular
   % velocities, but these ara actually not measured in rad/seg, because
   % the optical flow velocity isnt measured in pix/seg either!!! What it
   % gives us isjus ta difference. 
   %In which order do these rotations occur (to convert it to quaternions from euler)??? 
   % It looks like the Image velocity equations have already taken
   % the small angle asumption, and that is why it looks like the order
   % doesnt matter . Still... not sure about it. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%%  ------------ Test -
        
%--------------------2-D tests---------------------------------------------
%     OpticalFlow2D
%     Assume That the camera movement is uniform in the camera. Orientation is fixed. height is fixed. 

%     ODO_OpticalFlow2Dpos = [-mean(double(OpticalFlow_u(i,1:max_inliers(i))))*pixelsize*X(3)/f  mean(double(OpticalFlow_v(i,1:max_inliers(i))))*pixelsize*X(3)/f];
%     ODO_OpticalFlow2Dpos = [-mean(u)  -mean(v)];


%   tests movements in 2D
%     figure(1)
%     subplot(1,3,1)
%     title('Epipolar Search 5DoF');
%     actualpos = testpos_OF3D;
%     testpos_OF3D =  testpos_OF3D + ODO_OpticalFlowpos(1:2);
%     
%     plot([actualpos(1) testpos_OF3D(1) ],[actualpos(2) testpos_OF3D(2)])
%     hold on
%     % test remove outliers
%     
%     subplot(1,3,2)
%     title('Flow averaging 2DoF')
%     actualpos2 = testpos_OF2D;
%     testpos_OF2D =  testpos_OF2D + ODO_OpticalFlow2Dpos;
%     plot([actualpos2(1) testpos_OF2D(1) ],[actualpos2(2) testpos_OF2D(2)])
%     hold on
%     
%     subplot(1,3,3)
%     title('Optical Flow')
%     quiver(full_x, full_y, full_u, full_v, 'r')
%     hold on
%     quiver(x , y,  u , v, 'b')
%     xlabel('red = discarded outliers ; blue = inliers');
%     hold off
%     
%     subplot(2,2,4)
%     title('Histogram')
%     %u
%     [all_freq,all_value] = hist(full_u,20);
%     [inliers_freq,inliers_value] = hist(u,all_value);
%     stem(all_value,all_freq,'r');
%     hold on;
%     stem(inliers_value,inliers_freq,'b');
%     %v
%     [all_freq,all_value] = hist(full_v,20);
%     [inliers_freq,inliers_value] = hist(v,all_value);
%     stem(all_value,-all_freq,'r');
%     stem(inliers_value,-inliers_freq,'b');
%     hold off;
    
    
%     title('IMU integration')
%     actualpos3  = testpos_IMU;
%     testpos_IMU =  testpos_IMU + ODO_IMUpos(1:2);
%     % if abs(testpos_IMU - k(j,1:2)) > 1e-5
%     %     error('the do not coincide');
%     % else
%     %     fprintf('the IMU integrations coincide!! :):) \n');
%     % end
%     plot([actualpos3(1) testpos_IMU(1) ],[actualpos3(2) testpos_IMU(2)])
%     hold on
%     w = waitforbuttonpress;
% ---------------------------------------------------------------------

    
    
    
    
    
    

%% Kalman Filter
end


%-------------------------3-D Tests----------------------------------------

figure(5)
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
hold on
plot3(X(:,1),X(:,2),X(:,3),'r');
