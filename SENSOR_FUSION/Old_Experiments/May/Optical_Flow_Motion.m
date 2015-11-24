%Human and Pedestrian Nagivation homework
%Localization taskusing particle filters 
%Cristian Duguet 

% clc
close all
clear all

%Gravitational Constant
G =  9.80665;


%Careful. In this algorithm i did not take into account that the robot is
%rotating and hance, the angle depend on the orientation of the robot. This
%must be corrected in a future version.

experiment_date = '20120509_1920';


%% Parameters:
%in mm
sigmaIMUacc = [20 20 20];

sigmaIMUrot = [2 2 2];% in deg


%% Read files


ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));

load(strcat(experiment_date,'_Flow.txt.mat'));
load(strcat(experiment_date,'_PoseVicon.txt.mat'));


%We noticed that we have asynchronized video and IMU stream. So what we do
%is to make an algorithm which gets this data asynchronized??? TO SEE


%% IMU Construction
% Build an IMU sensor based on the Vicon position
 
IMUacc = diff(diff(PoseVicon(:,1:3)));
IMUacc = IMUacc + repmat(G*1000 0 0,size(IMUacc,1),1);



IMUrot = diff(PoseVicon(4:6));
IMUrot = IMUrot(1:end-1,:);

IMUTime = ViconTime(2:end-1);

% Build Noise




%% Optical Flow 

%utilize height Data:


% camera height: 









%Make the odometry data begin in (0,0,x)mm
 M{2} = M{2} - M{2}(1);
 M{3} = M{3} - M{3}(1);

 Time_instants = length(M{1});


 
 %ecuaciones de rigid motion 
 
 u = 1/Z *(-f *T_0 + x*T_2) + omega_0*(x*y/f) - omega_1 * (f + x^2/f) + omega_2*y;
 v = 1/Z *(-f *T_1 + y*T_2) + omega_0*(f + y^2/f) - omega_1 * (x*y/f) - omega_2*x;
 

%Odometry measurements
 % take the differential position for the odometer
 sigma_U = [ 30 30 30 ];
 %create odometry noise
 %we are using the normal distribution, but it can be any other    
    %2nd arg is position differences
    %3d argument is the repeated matrix of variances along K measurements
 U = random('Normal',[diff(M{2}) diff(M{3}) diff(M{4})],repmat([sigma_U],Time_instants-1,1));
 U = [[0 0 0] ; U];


%Measurement creation 
 Z = [0 0 ; diff(OF{2}) diff(OF{3})];
 % consider a variance for the measurements (mm):
 sigma_Z = [300 300 500];



%construct the actual measurements which wwill be evaluating the particles
RealZ = zeros(Time_instants,3);
RealZ(1,:) = [Z(1,1:2) 1000]; %mm         %initial height assumtion: 1m 

%Altimeter Construction
Alt = M{4};
sigma_Alt = [50];
 %create Altimeter noise
Alt = Alt + random('normal', 0,sqrt(sigma_Alt),size(Alt));


%Just U,sigma_U,Z1,Z2,Z3,Z4,sigma_Z,L1,L2,L3,L4 can be used to the position 
%estimation in the particle filter
%% Particle Filter Algorithm 
N_particles = 500;

%All particles are going to be in the matrix  (initialization)
X_P = repmat([0 0 1000],N_particles,1);

%as they are all on the same position, the first particles have the same
%weight
w_k = ones(N_particles,1)/N_particles;



% --------- Apartado para height estimation

% z: height
% f: focal length 
F = 6; %mm
pixelsize=3.6/480;
% -----------------------------------------------




% initialize posterior
POSTERIOR_MEAN = zeros(Time_instants,3);
POSTERIOR_MEAN(1,:)=X_P(1,:);
x_kminus2=POSTERIOR_MEAN(1,1);
y_kminus2=POSTERIOR_MEAN(1,2);

error = zeros(Time_instants,1);
Neff  = ones(Time_instants,1) * N_particles;




% Create Video
% aviobj = avifile(sprintf('OF-IMU-Alt_varZ_%d_varU_%d_varAlt_%d-%s',sigma_Z,sigma_U,sigma_Alt,datatoopen));
% aviobj.Quality = 100;

i=1;
hPlot= figure(1)
plot3(M{2}(1:i),M{3}(1:i),M{4}(1:i),'k','Linewidth',0.5);hold on;
plot3(POSTERIOR_MEAN(1:i,1),POSTERIOR_MEAN(1:i,2),POSTERIOR_MEAN(1:i,3),'b');
plot3(X_P(:,1),X_P(:,2),X_P(:,3),'r.');axis vis3d; axis square;
hold off
% aviobj = addframe(aviobj,hPlot);

display('Particle filter...GO!');
tic;

for i=2:Time_instants
%first iteration (draw x from proposal function) % it starts taking the
%odometry from U(2,:) and so on...


%% Sensor information 
%Note: do not confuse Z: sensor with z: position.

%IMPORTANT: Actual z position estimation: for estimating the height of the camera
%using the relationn between pixel movement (the Z sensor info) and actual movement, we use the
    %information las position and the one before that. We cannot use the
%odometry prior since that would make the information redundant. 
%The equation used is the pinhole equation
 
%IDEA: We could also use a larger memory to take the mean of z
%IDEA: We could also use separated particles for z

%Do not update if we do not move along x or y 

% RealZ(i,3) = POSTERIOR_MEAN(i-1,3) + 1/2* ((POSTERIOR_MEAN(i-1,1)-x_kminus2)*F/(Z(i-1,1) *pixelsize) + (POSTERIOR_MEAN(i-1,2)-y_kminus2)*F/(Z(i-1,2) *pixelsize));
% 
% if isnan(RealZ(i,3))
%     if Z(i-1,1)==0
%         if Z(i-1,2)==0
%             RealZ(i,3) = POSTERIOR_MEAN(i-1,3);
%         else 
%             RealZ(i,3) = POSTERIOR_MEAN(i-1,3) + (POSTERIOR_MEAN(i-1,2)-y_kminus2)*F/(Z(i-1,2) *pixelsize);
%         end
%     elseif Z(i-1,2)==0
%         RealZ(i,3) = POSTERIOR_MEAN(i-1,3) + (POSTERIOR_MEAN(i-1,1)-x_kminus2)*F/(Z(i-1,1) *pixelsize);
%     else
%         RealZ(i,3) = POSTERIOR_MEAN(i-1,3);
%     end
% end
% 
% RealZ(i,3) = 0.8* RealZ(i-1,3)+ 0.2*RealZ(i,3);

RealZ(i,3) = Alt(i); 

%Using this estimation of the height we calculate the movement using the
%pinhole equation too
  %this version uses the last posterior
RealZ(i,1) = POSTERIOR_MEAN(i-1,1) + RealZ(i,3) * ((Z(i,1)-1)*pixelsize)/F*2;
RealZ(i,2) = POSTERIOR_MEAN(i-1,2) + RealZ(i,3) * ((Z(i,2)+1)*pixelsize)/F*2;

  %this version use the last measurement position 
% RealZ(i,1) = RealZ(i-1,1)  + RealZ(i,3) * ((Z(i,1)-1)*pixelsize)/F*2;
% RealZ(i,2) = RealZ(i-1,2) + RealZ(i,3) * ((Z(i,2)+0.5)*pixelsize)/F*2;

%use them for the next iteration
x_kminus2 = POSTERIOR_MEAN(i-1,1);
y_kminus2 = POSTERIOR_MEAN(i-1,2);


%% Particle Updating

X_P = X_P + random('normal',repmat(U(i,:),N_particles,1),...
    repmat(sigma_U,N_particles,1));

%% Weight Update 
%We update the weights, evaluating the 2d multivariate normal pdfs, for n
%elements. We multiply the 4 results.


upd= mvnpdf(X_P, repmat([RealZ(i,1) RealZ(i,2) RealZ(i,3)],N_particles,1) ,...
    repmat(diag(sigma_Z),[1,1,N_particles])) ...
    .* mvnpdf(X_P(:,3), repmat(Alt(i),N_particles,1) ,...
    repmat(sigma_Alt,[1,1,N_particles]));




if upd==0
    fprintf('Zero weight on all particles...\n')
    fprintf('Trying bigger Z variance (4x)...\n');
    upd= mvnpdf(X_P, repmat([RealZ(i,1) RealZ(i,2) RealZ(i,3)],N_particles,1) ,...
        repmat(diag(4*sigma_Z),[1,1,N_particles]))...
        .* mvnpdf(X_P(:,3), repmat(Alt(i),N_particles,1) ,...
    repmat(sigma_Alt,[1,1,N_particles]));
    if upd==0
        fprintf('Zero weight on all particles...\n')
        fprintf('Trying bigger Z variance (10x)...\n');
        upd= mvnpdf(X_P, repmat([RealZ(i,1) RealZ(i,2) RealZ(i,3)],N_particles,1) ,...
            repmat(diag(10*sigma_Z),[1,1,N_particles]))...
            .* mvnpdf(X_P(:,3), repmat(Alt(i),N_particles,1) ,...
        repmat(sigma_Alt,[1,1,N_particles]));
        if upd==0
            fprintf('Incresing sensor variance does not seem to work.\n');
%              aviobj = close(aviobj);
            return;
        end
    end
end
        
    
w_k = w_k.* upd;
%normalize w_k
w_k = w_k/sum(w_k);

% if isnan(w_k)
% %         errordlg('The weights are NaN.','NaN');
%      cd(prev_fol); error('The weights are NaN.');
%     return;
% end

%Posterior
POSTERIOR_MEAN(i,:) = w_k'*X_P;



%% Resampling
%Effective number of particles
Neff(i) = 1/sum(w_k.^2);
if Neff(i) < 0.3*N_particles
    fprintf('Doing resampling on step %d...\n',i);
    %cumulative function of weights
    cum = cumsum(w_k);
    %uniformely distributed intervals
    u = linspace(0,1-1/N_particles,N_particles) + unifrnd(0, 1/N_particles);
    
    %initialize value asignation matrices
    X_P_new = zeros(N_particles,3);
    w_k_new = zeros(N_particles,1);
    for h = 1:N_particles
        index_found = find(cum>u(h),1);
        X_P_new(h,:) = X_P(index_found,:);
        w_k_new(h) = w_k(index_found);
    end
    X_P = X_P_new;
    w_k = w_k_new;
end




%plot figure in 3D
% hPlot=figure(1);
% plot3(M{2}(1:i),M{3}(1:i),M{4}(1:i),'k','Linewidth',0.5);hold on;
% plot3(POSTERIOR_MEAN(1:i,1),POSTERIOR_MEAN(1:i,2),POSTERIOR_MEAN(1:i,3),'r');
% plot3(X_P(:,1),X_P(:,2),X_P(:,3),'r.');
% plot3(RealZ(i,1),RealZ(i,2),RealZ(i,3),'bo','Linewidth',2);axis vis3d; axis equal;
% legend('Real Trajectory','Calculated trajectory','Particles','Measurement')
% hold off
% aviobj = addframe(aviobj,hPlot);

end
toc
%              aviobj = close(aviobj);
% display('Finished processing video.')

%% Plot Trajectories from Odometry and sensors:
hPlot=figure(1);
% plot real path
plot3(M{2}(1:i),M{3}(1:i),M{4}(1:i),'k','Linewidth',0.5);hold on;
% plot mean of the posterior + ESTIMATION 
plot3(POSTERIOR_MEAN(1:i,1),POSTERIOR_MEAN(1:i,2),POSTERIOR_MEAN(1:i,3),'r');
% plot particles
plot3(X_P(:,1),X_P(:,2),X_P(:,3),'r.');
% plot raw synthetic IMU odometry
URawPath = repmat([0 0 1000],Time_instants,1) + cumsum(U,1);
plot3(URawPath(:,1),URawPath(:,2),URawPath(:,3),'g');
% plot Sensor trajectory
plot3(RealZ(:,1),RealZ(:,2),RealZ(:,3),'b','Linewidth',2);axis vis3d; axis equal;
%legends and labels
legend('Real Path','Estimation','Particles','Odometry Raw','Optical Flow Raw')
xlabel('x'); ylabel('y')';zlabel('z');
hold off


figure(2)
subplot(1,2,1); 
plot(M{2},M{3}); 
subplot(1,2,2); plot(RealZ(:,1),RealZ(:,2));


figure(3)

plot(Neff); title('Evolution of effective number of particles'); grid on; 
% hold on;  plot(repmat(N_particles,1,K));
xlabel('time (steps)'); ylabel('Neff')

figure(4)
squareerror= sum((POSTERIOR_MEAN - [M{2} M{3} M{4}]).^2,2);
plot(squareerror);
% legend('x','y','z');
xlabel('time (steps)'); ylabel('squared error (mm2)')


%create annotations to plot
% sorted_w = sort(w_k,'descend');
% for j=1:N_particles
%     annotation('textbox',[ 0.5 0.6 0.1 0.05],...
%         'String', sprintf('%g',find(sorted_w == w_k(j))),...
%         'FontSize',8);
% end

