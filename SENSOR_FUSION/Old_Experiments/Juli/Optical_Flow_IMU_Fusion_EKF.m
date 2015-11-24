%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical Flow
% Cristian Duguet 

clc
close all
clear all

% Gravitational Constant
G =  9.80665;


% Experiment name 
experiment_date = '20120718_2021';
% experiment_date = '20120704_1913';

%% General Parameters:
%metrics in mm
sigma2IMUacc = [1 1 1]; % in (mm/s2)^2
sigma2IMUrotvel = deg2rad([5 5 5]).^2; % in rad^2


% Focal Length (mm):
f = 6; %mm 

% IMU Update Rate:
IMU_rate = 50; %Hz

% Scaling constant to be changed 
manual_scale = 36;

% Number of image instants which do not have Optical Flow
MAX_NOFLOW = 100;


%% Read files

%If the mat files are not created from the log files, then create them
try
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
catch ME
    ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
    ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));
    
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
end

%Analyze the frequency of messages  
delta = diff(ViconTime,1);
Viconfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

delta = diff(OpticalFlowTime);
OpticalFlowfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

% FIXME: sample rate of the camera does not seem to be 15FPS 
% display(sprintf('Mean rate of flor information: %g arrays per second \n',OpticalFlowfreq);


% Coverting ground truth into mm
ViconPose(:,1:3) = ViconPose(:,1:3)*1000;

% Pixel size 
pixel_size=3.75e-3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Preparation

%----------------------IMU Construction -----------------------------------

% Build an IMU sensor based on the Vicon Pose
timedifference = double(diff(ViconTime(:,1))) + double(diff(ViconTime(:,2)))/1e9; % sec
timedifference = repmat(timedifference,1,3);


% --- Rotation velocity: 
% The IMU is delivering the results in angular velocity


% Data Correction: We need to rotate the ViconPose by
% FIXME: Just for this experimen
%Rotation -90 in Z axis then 90 in X axis: [.5 -.5 -.5 .5]
ViconPose(:,4:7) = quatmult(repmat([.5 -.5 -.5 .5],size(ViconPose,1),1),ViconPose(:,4:7));
%FIXME: maybe its better to set the camera orientation as looking downwards
%([-1 0 0 0])

% Take the differential quaternion 
invPose = invquat(ViconPose(1:end,4:7)); %inverse quaternion
IMUrotvel = quatmult(ViconPose(2:end,4:7),invPose(1:end-1,:));
% IMUrotvel = quat_normalize(IMUrotvel);

%FIXME: Vicon sometimes deliver bad irregular rotations. Supress these
%making it equal to their previous rotation velocity
for i=1:size(IMUrotvel,1)
    if abs(IMUrotvel(i,4) - 1) > 0.1
        IMUrotvel(i,:) = [0 0 0 1];
        ViconPose(i+1,4:7) = ViconPose(i,4:7);
        invPose(i,:) = invPose(i-1,:);
        error('hola');
    end
end
%Take the angular velocity from the differential quaternion. 
IMUrotvel = 2*IMUrotvel(2:end,1:3)./timedifference(2:end,:);

% --- Acceleration: 
% Calculate velocity 

vel = diff(ViconPose(:,1:3)) ./timedifference; %mm/sec
IMUacc = diff(vel) ./timedifference(2:end,:) ; %mm/sec2

%%-- debug 
IMUacc_gt = IMUacc;
IMUacc_gt = quatmult([IMUacc_gt zeros(size(IMUacc_gt,1),1)],invPose(3:end,:));
IMUacc_gt = quatmult(ViconPose(3:end,4:7),IMUacc_gt);
IMUacc_gt = IMUacc_gt(:,1:3);
%%- debug


% Add Gravity
IMUacc = IMUacc + repmat([0 0 -G*1000],size(IMUacc,1),1);

% Rotate acceleration

IMUacc = quatmult([IMUacc zeros(size(IMUacc,1),1)],invPose(3:end,:));
IMUacc = quatmult(ViconPose(3:end,4:7),IMUacc);
IMUacc = IMUacc(:,1:3);


IMUTime_instants = size(IMUacc,1);

% ---- Noise: 
%FIXME; No bias considered
% Translational Noise 
% IMUacc = IMUacc + random('Normal',zeros(IMUTime_instants,3), repmat(sigma2IMUacc,IMUTime_instants,1));

% Rotational Noise
% IMUrotvel = IMUrotvel + random('Normal', zeros(IMUTime_instants,3), repmat(sigma2IMUrotvel,IMUTime_instants,1));

%-------------------------- Time adjustments ------------------------------

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

IMUTime = double((ViconTime(3:end,1))) + double((ViconTime(3:end,2)))/1e9 ;
OpticalFlowTime = double((OpticalFlowTime(:,1))) + double((OpticalFlowTime(:,2)))/1e9;

% the time in which both are transmiting
time_zero =  min(OpticalFlowTime(1), IMUTime(1));

OpticalFlowTime = double(OpticalFlowTime-time_zero);
IMUTime = double(IMUTime - time_zero);


% FIXME: there are still some errors in saving the files 
aux = max_feat;
max_feat = find(isnan(sum(OpticalFlow_u+OpticalFlow_v + OpticalFlow_x +OpticalFlow_y,1)),1)-1;
if isempty(max_feat), max_feat = aux;
end


OpticalFlow_u = OpticalFlow_u(:,1:max_feat);
OpticalFlow_v = OpticalFlow_v(:,1:max_feat);
OpticalFlow_x = OpticalFlow_x(:,1:max_feat);
OpticalFlow_y = OpticalFlow_y(:,1:max_feat);

for i = 1:size(matched_feat,1)
    
matched_feat(i) = min(matched_feat(i),max_feat);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Odometry and Predicion 

% ** Considerations:
% The imu data is going to be integrated according to the time steps of the
% Optical Flow (slower)
% In the firs iteration the IMU data will be integrated until its
% sample right after the new OpticalFlow sample timestamp. 
%  *IMPORTANT: While the imu is providing data in time differences (mm/s2),
%  the Optical Flow is considering velocity as just differences, not over a
%  time lapse. This means, if we want to integrate the IMU Data to the
%  Optical Flow, we have to consider  the rotational velocity as a Rotation
%  that has been made during a time slot, or a translation during that time
%  slot.

%% Counters and memory allocations
noflow_counter =0;

max_inliers = zeros(size(matched_feat));


%Sensors Information initialization
vel_IMUOF = zeros(length(OpticalFlowTime),3);
omega_IMU = zeros(length(OpticalFlowTime),3);

% % Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
% POSE_OF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison
% POSE_IMUOF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison


X_gt = repmat([ViconPose(1,1:3) 0 0 0 ViconPose(1,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1],length(OpticalFlowTime),1);
vel_gt = repmat( [0 0 0],length(OpticalFlowTime),1);
acc_gt = repmat( [0 0 0],length(OpticalFlowTime),1);

omega_m = zeros(length(OpticalFlowTime),3);
a_m = zeros(length(OpticalFlowTime),3);


%% --------------------------KALMAN FILTER---------------------------------

% we  express the first state as X0, the expected aas X_exp, the
% state itself as X 
% The error state as X_err, the prior as X_err_prior (no need to store in
% time)

%------------------------------Noises--------------------------------------

%Acceleration noise variance mm/s2
var_na = sigma2IMUacc;
sigma2IMUrotvel = deg2rad([5 5 5]).^2; % in rad^2
% Acceleration Bias noise variance mm/s2
sigma_nba = [0.5 0.5 0.5];        var_nba = sigma_nba.^2;

% Rotational velocity noise variance
var_nomega = sigma2IMUrotvel;
%Rotational velocity bias noise variance
sigma_nbomega = [0.1 0.1 0.1];    var_nbomega = sigma_nbomega.^2;

%scale noise
var_nL = 0.2^2; % relative to unit

% Covariance Matrix of the Optical Flow measurements
R = diag([5 5 5].^2); %mm2 
%FIXME: better representation

% Noises needed by Stephan weiss q_d calculationa lgorithm 
    var_nqvw = [0 0 0];
    var_nqci = [0 0 0];
    var_npic = [0 0 0];


%-----------------------State Initialization-------------------------------

%We assume the initial position is at one meter height
p_w_i   = ViconPose(3,1:3);   %FIXME: We are taking the ground truth initial position, for comparison
v_w_i   = vel(2,:); % Assume zero velocity  %FIXME
q_w_i   = [0 0 0 1];  % Assume no rotation at the beginning 
b_omega = [0 0 0];
b_a     = [0 0 0];
L       = 1;
p_i_c   = [0 -4 -4];   %Calibration state: Distance from IMU to camera
q_i_c   = [0 0 0 1];     %calibration state: Orientation of camera w.r.t IMU frame 
 


% For the initial state, the state vector wil be 
X0 = [ p_w_i v_w_i q_w_i b_omega b_a L p_i_c q_i_c] ;  % 24 element vector 
%Expected state vector 

%State vector itself
X = repmat(X0,length(OpticalFlowTime),1);
X_exp = repmat(X0,length(OpticalFlowTime),1);

% Error state vector
Deltap_w_i = [0 0 0];
Deltav_w_i = [0 0 0];
deltaomega_w_i = [0 0 0];
Deltab_w   = [0 0 0];
Deltab_a   = [0 0 0];
DeltaL     = 0;
Deltap_i_c = [0 0 0];
deltaomega_i_c = [0 0 0];

%And the error state 
X_err0 = zeros(1,22);
X_err = repmat(X_err0,length(OpticalFlowTime),1);

% For the initial state, the error state covariance matrix will be
% trueerror = zeros(1,22);
% trueerror(1:6) = X0(1:6)-X_gt(1:6);
% aux = quatmult(X0(7:10),invquat(X_gt(6:10)));
% trueerror(7:9) = aux(1:3);
% trueerror(10:19)=X0(11:20)-X_gt(11:20);
% aux = quatmult(X0(21:24),invquat(X_gt(21:24)));
% trueerror(20:22) = aux(1:3);
% P_kminusone_kminusone = trueerror'*trueerror +.5*ones(22);
% P_kminusone_kminusone = eye(22);

P_kminusone_kminusone = [0.016580786012789, 0.012199934386656, -0.001458808893504, 0.021111179657363, 0.007427567799788, 0.000037801439852, 0.001171469788518, -0.001169015812942, 0.000103349776558, -0.000003813309102, 0.000015542937454, -0.000004252270155, -0.000344432741256, -0.000188322508425, -0.000003798930056, 0.002878474013131, 0.000479648737527, 0.000160244196007, 0.000012449379372, -0.000025211583296, -0.000029240408089, -0.000001069329869, -0.001271299967766, -0.000133670678392, -0.003059838896447
		, 0.012906597122666, 0.050841902184280, -0.001973897835999, 0.017928487134657, 0.043154792703685, 0.000622902345606, 0.002031938336114, 0.000401913571459, -0.000231214341523, -0.000016591523613, 0.000011431341737, 0.000007932426867, 0.000311267088246, -0.000201092426841, 0.000004838759439, 0.008371265702599, -0.000186683528686, 0.000139783403254, 0.000070116051011, -0.000021128179249, -0.000028597234778, -0.000006006222525, -0.002966959059502, 0.000313165520973, 0.003179854597069
		, -0.001345477564898, -0.000886479514041, 0.014171550800995, -0.002720150074738, 0.005673098074032, 0.007935105430084, 0.000687618072508, 0.000684952051662, 0.000022000355078, -0.000008608300759, -0.000000799656033, 0.000001107610267, -0.000106383032603, -0.000356814673233, -0.000068763009837, -0.000051146093497, -0.000091362447823, 0.000293945574578, -0.000256092019589, 0.000042269002771, -0.000009567778418, -0.000017167287470, 0.004592386869817, -0.001581055638926, 0.000227387610329
		, 0.020963436713918, 0.016241565921214, -0.002606622877434, 0.043695944809847, 0.008282523689966, -0.001656117837207, 0.001638402584126, -0.002060006975745, -0.001362992588971, -0.000001331527123, 0.000032032914797, 0.000004134961242, 0.000341541553429, -0.000100600014193, 0.000025055557965, 0.003723777310732, -0.000161259841873, 0.000175908029926, -0.000010843973378, -0.000001022919132, -0.000020982262562, -0.000009716850289, -0.002231080476166, -0.001033766890345, -0.003628168927273
		, 0.009314922877817, 0.046059780658109, 0.003565024589881, 0.015262116382857, 0.065035219304194, -0.001635353752413, 0.002492076189539, 0.001255538625264, -0.000034886338628, -0.000029672138211, 0.000006695719137, 0.000006779584634, 0.000273857318856, 0.000241559075524, 0.000026819562998, 0.007341077421410, -0.000245364703147, -0.000214640089519, 0.000072765069578, -0.000031941424035, 0.000014164172022, -0.000014177340183, -0.000530959567309, 0.000080230949640, 0.003376885297505
		, -0.000029025742686, 0.000535037190485, 0.007958782884182, -0.001871298319530, -0.002083832757411, 0.012983170487598, 0.000132746916981, 0.000083483650298, 0.000020140288935, -0.000001280987614, 0.000000838029756, -0.000000023238638, -0.000309256650920, 0.000094250769772, -0.000143135502707, 0.000262797080980, 0.000133734202454, 0.000025809353285, 0.000051787574678, 0.000002954414724, -0.000012648552708, -0.000004097271489, 0.002381975267107, -0.001036906319084, 0.000115868771739
		, 0.001237915701080, 0.002441754382058, 0.000642141528976, 0.001714303831639, 0.003652445463202, 0.000133021899909, 0.000491964329936, 0.000029132708361, 0.000054571029310, -0.000003531797659, 0.000002108308557, -0.000000655503604, -0.000036221301269, -0.000080404390258, -0.000002011184920, 0.000409618760249, 0.000006455600111, 0.000037893047554, 0.000004332215700, -0.000003727533693, 0.000000308858737, -0.000004128771100, 0.000121407327690, -0.000116077155506, -0.000044599164311
		, -0.001129210933568, 0.000810737713225, 0.000687013243217, -0.002320565048774, 0.001923423915051, 0.000083505758388, 0.000045906211371, 0.000464144924949, -0.000074174151652, -0.000001593433385, -0.000002820148135, 0.000001999456261, 0.000068256370057, -0.000050158974131, -0.000000228078959, 0.000046796063511, -0.000043197112362, 0.000007902785285, 0.000000020609692, 0.000001805172252, 0.000002146994103, 0.000005750401157, 0.000309103513087, 0.000176510147723, 0.000423690330719
		, 0.000118011626188, -0.000151939328593, -0.000003895302246, -0.001370909458095, 0.000050912424428, 0.000014452281684, 0.000048567151385, -0.000077773340951, 0.000550829253488, -0.000001499983629, -0.000001785224358, -0.000005364537487, 0.000036601273545, 0.000003384325422, -0.000000535444414, -0.000032994187143, -0.000004973649389, -0.000005428744590, 0.000002850997192, -0.000006378420798, -0.000000001181394, -0.000014301726522, 0.000038455607205, 0.000110350938971, -0.000142528866262
		, -0.000005270401860, -0.000021814853820, -0.000010366987197, -0.000002004330853, -0.000038399333509, -0.000001674413901, -0.000004404646641, -0.000002139516677, -0.000001756665835, 0.000002030485308, -0.000000003944807, 0.000000005740984, 0.000000210906625, 0.000000302650227, 0.000000014520529, -0.000003266286919, -0.000000158321546, -0.000000508006293, -0.000000000135721, -0.000000498539464, 0.000000163904942, 0.000000129053161, -0.000003222034988, 0.000000064481380, -0.000001109329693
		, 0.000016356223202, 0.000012074093112, -0.000001861055809, 0.000034349032581, 0.000006058258467, 0.000000706161071, 0.000001988651054, -0.000003017460220, -0.000001874017262, -0.000000012182671, 0.000002030455681, -0.000000019800818, 0.000000488355222, 0.000001489016879, 0.000000028100385, 0.000002786101595, -0.000000046249993, 0.000000097139883, 0.000000389735880, -0.000000195417410, -0.000000460262829, 0.000000210319469, -0.000002235134510, -0.000002851445699, -0.000002883729469
		, -0.000003154072126, 0.000010432789869, 0.000002047297121, 0.000005626984656, 0.000009913025254, 0.000000398401049, -0.000000326490919, 0.000002058769308, -0.000005291111547, 0.000000001086789, 0.000000001772501, 0.000002006545689, 0.000000044716134, 0.000000414518295, -0.000000135444520, 0.000001531318739, -0.000000211673436, 0.000000405677050, -0.000000796855836, -0.000000266538355, -0.000000133632439, -0.000000338622240, -0.000000150597295, -0.000000563086699, 0.000003088758497
		, -0.000348907202366, 0.000314489658858, -0.000097981489533, 0.000332751125893, 0.000276947396796, -0.000311267592250, -0.000035302086269, 0.000070545012901, 0.000036626247889, 0.000000400828580, 0.000000087733422, 0.000000120709451, 0.001026573886639, 0.000013867120528, 0.000031828760993, 0.000009746783802, -0.000458840039830, -0.000019468671558, -0.000043520866307, 0.000007245947338, 0.000003901799711, -0.000004201599512, -0.000047176373840, 0.000119567188660, 0.000003684858444
		, -0.000190283000907, -0.000192352300127, -0.000359131551235, -0.000107453347870, 0.000258576553615, 0.000091496162086, -0.000081280254994, -0.000048304910474, 0.000002800928601, 0.000000908905402, 0.000001125333299, 0.000000471832044, 0.000019874619416, 0.001029579153516, 0.000011053406779, 0.000021449316681, 0.000016006639334, -0.000412772225495, 0.000006993477540, 0.000002648721730, 0.000004792699830, -0.000004141354722, -0.000083992422256, 0.000015935718681, -0.000000338251253
		, -0.000004368584055, 0.000003124910665, -0.000067807653083, 0.000024474336501, 0.000022105549875, -0.000144033820704, -0.000002164571960, -0.000000083713348, -0.000000674226005, 0.000000019237635, 0.000000025526504, -0.000000057252892, 0.000032366581999, 0.000010736184803, 0.000111095066893, 0.000000615680626, -0.000015341510438, -0.000007700695237, -0.000023026256094, 0.000000638926195, 0.000000960343604, 0.000000817586113, -0.000026575050709, 0.000013993827719, -0.000002316938385
		, 0.002973222331656, 0.008292388147295, -0.000211655385599, 0.003951267473552, 0.006718811356807, 0.000277369882917, 0.000349425829596, -0.000014812000602, -0.000045952715508, -0.000002513020002, 0.000002692914948, 0.000001078825296, 0.000009897987444, 0.000020034595279, 0.000000809851157, 0.001554211174363, 0.000023959770856, -0.000037670361809, -0.000009320812655, -0.000004598853139, -0.000006284196194, -0.000000693801636, -0.000469324632849, 0.000014818785588, 0.000277219840791
		, 0.000476557664133, -0.000191539372645, -0.000089666716294, -0.000163721235917, -0.000235017605089, 0.000134712473215, 0.000007671308678, -0.000041648250772, -0.000005375975547, 0.000000156986772, 0.000000504340505, -0.000000198574002, -0.000458130878121, 0.000014584188938, -0.000015616513739, 0.000023678958593, 0.000535136781135, -0.000016449781236, 0.000040831795426, -0.000013702650244, -0.000000627377616, -0.000004196881223, 0.000002230529685, -0.000050724631819, -0.000004714535751
		, 0.000162219848991, 0.000116427796874, 0.000292562152669, 0.000173404902614, -0.000249216364740, 0.000026816594889, 0.000036367682776, 0.000005763510102, -0.000005320926337, -0.000000071291000, -0.000000112152457, 0.000000334342568, -0.000022684595881, -0.000410859858969, -0.000007890929454, -0.000040454023111, -0.000011131820455, 0.000458907544194, -0.000005285694195, 0.000002246982110, -0.000002222041169, 0.000001951461640, 0.000047488638766, -0.000029510929794, 0.000005816436594
		, 0.000010794825884, 0.000058045653749, -0.000260506684499, -0.000007544850373, 0.000048451414581, 0.000048500128303, 0.000002555777025, -0.000001118968589, 0.000001725856751, 0.000000113523451, 0.000000356160739, -0.000000287211392, -0.000041197824317, 0.000004749859562, -0.000021745597805, -0.000011794173035, 0.000040317421040, -0.000001104681255, 0.000325476240984, 0.000006084247746, -0.000006253095726, -0.000005627495374, 0.000013663440542, -0.000012536337446, 0.000000477230568
		, -0.000028222744852, -0.000029726624789, 0.000042365440829, -0.000004529013669, -0.000041974513687, 0.000002547416367, -0.000004149622895, 0.000001656132079, -0.000006464228083, -0.000000593440587, -0.000000063566120, -0.000000230872869, 0.000007212338790, 0.000002222629345, 0.000000642817161, -0.000006111733946, -0.000013813495990, 0.000002643879751, 0.000005887006479, 0.000020142991502, -0.000000692093175, -0.000000188761575, 0.000017519903352, -0.000002456326732, 0.000001576856355
		, -0.000026132063406, -0.000024675067133, -0.000008452766004, -0.000014350608058, 0.000014404004024, -0.000011620075371, 0.000000539065468, 0.000001829895964, -0.000000462890431, 0.000000223093202, -0.000000499925964, -0.000000094710754, 0.000003954308159, 0.000004249241909, 0.000000876422290, -0.000005419924437, -0.000001021458192, -0.000002052781175, -0.000007397128908, -0.000000347703730, 0.000021540076832, 0.000001455562847, 0.000005351749933, 0.000020079632692, 0.000006997090317
		, 0.000001606076924, 0.000001031428045, -0.000015843471685, -0.000005357648114, -0.000007152430254, -0.000003359339850, -0.000003466742259, 0.000005980188844, -0.000014512044407, 0.000000136766387, 0.000000188396487, -0.000000299050190, -0.000004280062694, -0.000005018186182, 0.000000751147421, 0.000000382366121, -0.000004319412270, 0.000002858658354, -0.000005774838189, -0.000000199234914, 0.000001477444848, 0.000021955531390, -0.000005912741153, 0.000006848954650, 0.000000718992109
		, -0.001250410021685, -0.002465752118803, 0.004640769479530, -0.002397333962665, 0.000543954908379, 0.002370095810071, 0.000159513911164, 0.000327435894035, 0.000051354259180, -0.000002658607585, -0.000001766738193, -0.000000182288648, -0.000049404478395, -0.000084546262756, -0.000026628375388, -0.000398670523051, 0.000000139079122, 0.000048715190023, 0.000014902392001, 0.000017378375266, 0.000005675773452, -0.000005943594846, 0.013030218966816, 0.002362333360404, 0.000426396397327
		, -0.000130856879780, 0.000387010914370, -0.001570485481930, -0.001207751008090, 0.000021063199750, -0.001030927710933, -0.000109925957135, 0.000181001368406, 0.000107869867108, 0.000000177851848, -0.000002935702240, -0.000000493441232, 0.000119019560571, 0.000014103264454, 0.000013824858652, 0.000027253599949, -0.000051452899775, -0.000028435304764, -0.000013422029969, -0.000002043413021, 0.000020290127027, 0.000006914337519, 0.002362694187196, 0.016561843614191, 0.000974154946980
		, -0.002974278550351, 0.003344054784873, 0.000125156378167, -0.003468124255435, 0.003442635413150, 0.000109148337164, -0.000076026755915, 0.000385370025486, -0.000148952839125, -0.000000760036995, -0.000002603545684, 0.000003064524894, 0.000001812974918, -0.000002381321630, -0.000002469614200, 0.000309057481545, -0.000004492645187, 0.000007689077401, 0.000001009062356, 0.000001877737433, 0.000007317725714, 0.000000467906597, 0.000372138697091, 0.000966188804360, 0.011550623767300;];
    P_kminusone_kminusone = P_kminusone_kminusone(1:22,1:22);

%FIXME




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First iteration

% We process first the IMU information. The received Optical Flow data is 
% ignored until the first IMU data has arrived. 

OpticalFlowstart = find(OpticalFlowTime > IMUTime(1),1);
IMUsecondit = find(IMUTime > OpticalFlowTime(OpticalFlowstart),1);

delta_t = IMUTime(IMUsecondit) - IMUTime(1);


% ----------------------IMU Strapdown algorithm----------------------------
 
% The initial velocity and position information is in X0

[X_exp(1,:),omega_m(1,:),a_m(1,:)] = IMU_Strapdown(X0,IMUrotvel(1:IMUsecondit,:),IMUacc(1:IMUsecondit,:),IMUTime(1:IMUsecondit));

% -------------------------- Optical Flow----------------------------------
% Nothing to do here

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

%% %%%%%%%%%%%%%%%%%%%%%%%Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------- Expected State Values -----------------------------

p_w_i_hat       = X_exp(1,1:3);
v_w_i_hat       = X_exp(1,4:6);
q_w_i_hat       = X_exp(1,7:10);  % Expectation fof the quaternion from world to IMU 
% q_w_i_hat = quatmult(ODO_IMUrot,X(7:10)); 
% q_w_i_hat = quat_normalize(q_w_i_hat);
b_omega_hat     = X_exp(1,11:13); % Expectation of the bias of the IMU gyro measurements
b_a_hat         = X_exp(1,14:16); % Expectation of the bias of the IMU accelerometer measurements
L_hat           = X_exp(1,17);
p_i_c_hat       = X_exp(1,18:20);
q_i_c_hat       = X_exp(1,21:24);
% FIXME: b_a_hat, b_omega_hat, p_i_c_hat and q_i_c_hat are not being changed in the
% expectancy, with respect to the state K

% --------------------Construct variables from expected values ----------

% Cross product matrices from the unbiased measurements 
Cross_a_x = VectorToCrossMat(a_m(1,:) - b_a_hat);
Cross_omega_x = VectorToCrossMat(omega_m(1,:) - b_omega_hat);
Cross_omega_x2 = Cross_omega_x*Cross_omega_x;

%Make rotation matrix from the quaternion 
C_q_w_i_hat = QuatToRotMat(q_w_i_hat);
C_q_i_c_hat = QuatToRotMat(q_i_c_hat);

% Make cross product matrices from v_w_i_hat and p_i_c_hat
Cross_v_w_i_hat = VectorToCrossMat(v_w_i_hat);
Cross_p_i_c_hat = VectorToCrossMat(p_i_c_hat);

%--------------------------------------------------------------------------

% ---------------------------State Matrices--------------------------------

F_d  = ConstructFd(C_q_w_i_hat,Cross_a_x,Cross_omega_x,Cross_omega_x2,delta_t);

%%------------------- Qd Calculation: Numerical -------------------------

% [F_d,Q_d_num] = GetMatricesNumerical(C_q_w_i_hat, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)

%%-------------------------------------------------------------------------

%%------------------- Qd Calculation: Analytical --------------------------

% Q_d= ConstructQdAnalytic(delta_t,C_q_w_i_hat, Cross_a_x,Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL);
%%-------------------------------------------------------------------------

%%-----------------------Stephan Weiss approach ----------------------------
Qd_weiss = calc_Q(delta_t,q_w_i_hat,omega_m(1,:) - b_omega_hat,a_m(1,:) - b_a_hat,...
    var_na,var_nba,...
    var_nomega,var_nbomega,...
    var_nL,...
    var_nqvw,...
    var_nqci,...
    var_npic);

Q_d = Qd_weiss(1:22,1:22);
%%-------------------------------------------------------------------------
%-----------------------PROPAGATION/PREDICTION-----------------------------

P_k_kminusone = F_d*P_kminusone_kminusone*F_d' + Q_d;
%bibliography Maybeck cap 4.9


%%-------------------------UPDATE/MEASUREMENT------------------------------

% It is important to notice that we just get measurements in the velocity.
% Then, we only write the H matrix for velocity measurements 


% Construct Measurement Matrix for velocity 

H = ConstructHv(L_hat,C_q_i_c_hat,C_q_w_i_hat,Cross_v_w_i_hat,Cross_omega_x,Cross_p_i_c_hat,p_i_c_hat,v_w_i_hat);

% Compute the residual in velocity 
r = [0 0 0]';  % For this first iteration we do not have Visual measurements 

% Covariance Matrix of the Optical Flow measurements
R; % defined in initialization of Kalman Filter

%Compute the innovation 
S = H * P_k_kminusone * H' + R;

% Compute Kalman Gain 
K  = P_k_kminusone*H'/S;

% Compute correction
X_correction = K*r;

% Update Coviariance Matrix 
P_k_k =  (eye(22) - K*H) * P_k_kminusone * (eye(22) - K*H)' + K*R*K';

% Update Error State Vector
X_err(1,:) = X_correction'; % again, the state is expressed horizontal

% Update State Vector 
X(OpticalFlowstart,1:6) = X_exp(1,1:6) + X_err(1,1:6);
X(OpticalFlowstart,7:10) = quatmult(quat_normalize([0.5*X_err(1,7:9) 1]), X_exp(1,7:10));
X(OpticalFlowstart,11:20) = X_exp(1,11:20) + X_err(1,10:19);
X(OpticalFlowstart,21:24) = quatmult(quat_normalize([0.5*X_err(1,20:22) 1]), X_exp(1,21:24));


X_err(OpticalFlowstart,:) = X_err(1,:);
%FIXME: we are not using the states in time instants 2:OpticalFlowTime-1

X_exp(OpticalFlowstart,:) = X_exp(1,:);





%%-debug
%FIXME: this is a debug: Ground truth
acc_gt(1,1:3) = IMUacc_gt(IMUsecondit,:);
vel_gt(1,1:3) = vel(IMUsecondit+1,:);
X_gt(1,:) = [ViconPose(IMUsecondit+2,1:3) vel_gt(1,1:3) ViconPose(IMUsecondit+2,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1];

fprintf('Result of the Kalman Filter: \n');
fprintf('%d:X0_vel \t X_exp_vel  \t Measurement \t Meas_exp \t X_corr_vel \t X_post_vel \t X_gt_vel \n',1);
for ii=1:3
    fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
             X0(3+ii) ,X_exp(1,3+ii), 0    ,0   ,  X_err(OpticalFlowstart,3+ii)  ,    X(OpticalFlowstart,3+ii) ,  X_gt(1,3+ii) );    
end
%%-debug

% TODO: use IMU acceleration as a orientation measurement 

%%-------------------------------------------------------------------------
%% Second iteration and So on: here we include the Optical Flow

index_end = IMUsecondit;

for i = OpticalFlowstart+1: length(OpticalFlowTime)
    
    %Time interval
    delta_t = OpticalFlowTime(i) - OpticalFlowTime(i-1);
    
    %% ================= IMU Strapdown Algorithm ==========================
    index_start = index_end; %use the end index of last iteration
    index_end = find(IMUTime > OpticalFlowTime(i),1); % FIXME: can crash if theres no IMU sample after the last OF sample
    
    
    % The initial velocity and position information is in X(i-1,:)
    %IMPORTANT: We are including the elements index_start(actual) =
    %index_end(previous) because wee are considering the timebetween two
    %iterations. This means, we just use IMUTime(indexstart). Anyways, it
    %is done so because of symmetryp
%     X(i-1,7:10) = ViconPose(index_start,4:7);
    [X_exp(i,:),omega_m(i,:),a_m(i,:)] = IMU_Strapdown(X(i-1,:),IMUrotvel(index_start:index_end,:),IMUacc(index_start:index_end,:),IMUTime(index_start:index_end));

    %% =====================OPTICAL FLOW===================================
    ama=1/36; % FIXME : get rid of this constant
    
    % Set the pixel position referenced to the center
    x = (OpticalFlow_x(i,1:matched_feat(i)))'-Image_width/2+0.5;
    y = (OpticalFlow_y(i,1:matched_feat(i)))'-Image_height/2+0.5;
    x = x *(ama/f); y = -y *(ama/f);
    z = ones(size(x));
    
    u = (OpticalFlow_u(i,1:matched_feat(i)))'*(ama/f);
    v = -(OpticalFlow_v(i,1:matched_feat(i)))'*(ama/f);
    
    %FIXME now they are just being divided by f, but not multiplied by
    %pixelsize
    
    
    %     %%In case of no flow detected
    %     if length(u) < 5
    %         if     noflow_counter > MAX_NOFLOW
    %             error('Optical Flow information too poor. Need more features.\n');
    %             return
    %         end
    %         warning('No flow detected in this image frame.\n');
    %         noflow_counter  = noflow_counter+1;
    %         %keep constant movement
    %         POSE_OF(i,:) = [ POSE_OF(i-1,1:3)+ODO_OpticalFlowpos quatmult(ODO_OpticalFlowrot,POSE_OF(i-1,4:7))];
    %         if length(u) < 2
    %             POSE_IMUOF = [ POSE_IMUOF(i-1,1:3)+ODO_OpticalFlowpos quatmult(ODO_OpticalFlowrot,POSE_IMUOF(i-1,4:7))];
    %         end
    %     continue
    %     end

    
    %% ---------------------RANSAC-----------------------------------------
    
    
    
    % First method: Calculating the Homogrpahy Matrix from a random subset
    % of sample_num points of all the matched features in the optical flow.
    % After max_Iter iterations, return the inliers of the most accepted
    % model (The model which can include the most inliers).
    % max_Iter = 1000;
    % sample_num = 5;
    % RANSAC_Homography_treshold = 1;
    %
    % inliers_idx = RANSAC_Homography(x,y,u,v,max_Iter,sample_num,RANSAC_Homography_treshold);
    
    
    %----------------------Linear RANSAC-----------------------------------
    % Second method: Matching the Consensus of the random samples in a
    % linear model
    
    max_Iter = 1000;
    sample_num = 50;
    RANSAC_Linear_treshold = 1;
    
    inliers_idx = RANSAC_Linear(u,v,max_Iter,sample_num,RANSAC_Linear_treshold);
    
    if length(inliers_idx) < 5
        warning('Inliers information too poor. Avoiding RANSAC selection.\n');
        w = zeros(size(u));
    else
        full_x=x;
        full_y=y;
        full_u=u;
        full_v=v;
        x = x(inliers_idx);
        y = y(inliers_idx);
        z = z(inliers_idx);
        u = u(inliers_idx);
        v = v(inliers_idx);
        w = zeros(size(u));
        
        max_inliers(i) = length(inliers_idx);
    end
    
    %%----------------RANSAC_END-------------------------------------------
    
    %% ---------------EPIPOLAR RECONSTRUCTION- 5DoF-------------------------
    %     Case 5DoF
    %     [omega5DoF(i,:),vel5DoF(i,:),correl5DoF] = Continuous8point5DoF(x,y,z,u,v,w);
    
    % Case 3DoF w=(0,0,0) FIRST UNROTATE!!!
    % As we get the rotation information from the IMU, we restrict the
    % solution of the OF epipolar search to the search of the translational
    % velocity
    
    %% ------------ SENSOR COLLIGATION and EPIPOLAR RECONSTRUCTION --------
    % Unrotate Optical Flow for use in 2DoF Epipolar Algorithm
    % The rotation velocity of the IMU has been integrated over time. The
    % result of it is the rotation quaternion. As the Optical Flow here is
    % a distance vector rather than a velocity vector, we use the
    % integrated rotation of the IMU to unrotate the Optical Flow
    
    rotation  = omega_m(i,:) * (IMUTime(index_end) - IMUTime(index_start+1));
    
    % Assuming smallangle:
    unrotation_OF = cross(repmat(rotation,max_inliers(i),1),[x y z],2);
    
    %UNROTATE
    u = u + unrotation_OF(:,1);
    v = v + unrotation_OF(:,2);
    w = w + unrotation_OF(:,3);
    
    %Calculate using the Epipolar Constraint over the Optical Flow, and the
    %Scale state
    [vel_IMUOF(i,:)] = 100*Continuous8point2DoF(x,y,z,u,v,w);
%     [vel_IMUOF(i,:)] = mean(vel(index_start+2:index_end+1,:),1);
    
    
    %% --------------------------SCALE RECOVERY -------------------------------
    % NOT IN USE: SCALE RESOLVED IN SENSOR FUSION FILTER
    
    %      crossprod = cross(repmat(ODO_IMUrot,1,max_inliers(i)), [x';y';z'],1);
    %Here we consider that the Optical Flow is unrotated
    %      M = [[diag(u), diag(x), -repmat(vel0(1),max_inliers(i),1)];
    %           [diag(v), diag(y), -repmat(vel0(2),max_inliers(i),1)];
    %           [diag(w), diag(z), -repmat(vel0(3),max_inliers(i),1)]];
    %
    %         [~,~,V2] = svd(M);
    %         lambda_scale = V2(:,end);
    %         manual_scale = lambda_scale(end);
    
    %FIXME: Using the unifying scale right??
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %--------------------2-D tests---------------------------------------------
    %     OpticalFlow2D
    %     Assume That the camera movement is uniform in the camera. Orientation is fixed. height is fixed.
    
    %     ODO_OpticalFlow2Dpos = [-mean(double(OpticalFlow_u(i,1:max_inliers(i))))*pixelsize*POSE(3)/f  mean(double(OpticalFlow_v(i,1:max_inliers(i))))*pixelsize*POSE(3)/f];
    %     ODO_OpticalFlow2Dpos = [-mean(u)  -mean(v)];
    
    
    %   tests movements in 2D
    %     figure(1)
    %     subplot(1,3,1)
    %     title('Epipolar Search 5DoF');
    %     actualpos = testpos_OF3D;
    %     testpos_OF3D =  testpos_OF3D + ODO_OpticalFlowpos(1:2);
    %
    %     plot([actualpos(1) testpos_OF3D(1) ],[actualpos(2) testpos_OF3D(2)])
    %     hold on[zeros(7,13)]];
    
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
    %     xlabel('red = discarded outliers ; blue[zeros(7,13)]];
    
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
    %     [inliers_freq,inliers_value] = hist(v,[zeros(7,13)]];
    
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
    %     hold on[zeros(7,13)]];
    
    %     w = waitforbuttonpress;
    % ---------------------------------------------------------------------
        
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Error State Covariance Matrix 
    P_kminusone_kminusone = P_k_k;    

    % --------------------- Expected State Values -----------------------------
    
    p_w_i_hat       = X_exp(i,1:3);
    v_w_i_hat       = X_exp(i,4:6);
    q_w_i_hat       = X_exp(i,7:10);  % Expectation of the quaternion from world to IMU
    % q_w_i_hat = quatmult(ODO_IMUrot,X(7:10));
    % q_w_i_hat = quat_normalize(q_w_i_hat);
    b_omega_hat     = X_exp(i,11:13); % Expectation of the bias of the IMU gyro measurements
    b_a_hat         = X_exp(i,14:16); % Expectation of the bias of the IMU accelerometer measurements
    L_hat           = X_exp(i,17);
    p_i_c_hat       = X_exp(i,18:20);
    q_i_c_hat       = X_exp(i,21:24);
    % FIXME: b_a_hat, b_omega_hat, p_i_c_hat and q_i_c_hat are not being changed in the
    % expectancy, with respect to the state K
    
    
    % --------------------Construct variables from expected values ----------
    
    
    
    % Cross product matrices from the unbiased measurements
    Cross_a_x = VectorToCrossMat(a_m(i,:) - b_a_hat);
    Cross_omega_x = VectorToCrossMat(omega_m(i,:) - b_omega_hat);
    Cross_omega_x2 = Cross_omega_x*Cross_omega_x;
    
    %Make rotation matrix from the quaternion
    C_q_w_i_hat = QuatToRotMat(q_w_i_hat);
    C_q_i_c_hat = QuatToRotMat(q_i_c_hat);
    
    % Make cross product matrices from v_w_i_hat and p_i_c_hat
    Cross_v_w_i_hat = VectorToCrossMat(v_w_i_hat);
    Cross_p_i_c_hat = VectorToCrossMat(p_i_c_hat);
    
    %--------------------------------------------------------------------------
    
    
    % ---------------------------State Matrices--------------------------------
    
    F_d  = ConstructFd(C_q_w_i_hat,Cross_a_x,Cross_omega_x,Cross_omega_x2,delta_t);
    
    %%------------------- Qd Calculation: Numerical -------------------------
    
%     [Q_d_num] = GetMatricesNumerical(C_q_w_i_hat, Cross_a_x, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL, delta_t);
    
    %%-------------------------------------------------------------------------
    
    %%------------------- Qd Calculation: Analytical --------------------------
    
%   Qd= ConstructQdAnalytic(delta_t,C_q_w_i_hat, Cross_omega_x,Cross_omega_x2, var_na, var_nba, var_nomega, var_nbomega, var_nL)

    %%-------------------------------------------------------------------------
    
    %%-----------------------Stephan Weiss approach ----------------------------
    
    Qd_weiss = calc_Q(delta_t,q_w_i_hat,omega_m(i,:) - b_omega_hat,a_m(i,:) - b_a_hat,...
        var_na,var_nba,...
        var_nomega,var_nbomega,...
        var_nL,...
        var_nqvw,...
        var_nqci,...
        var_npic);
    
    Q_d = Qd_weiss(1:22,1:22);
%     
    %%-------------------------------------------------------------------------
    
    %%-------------------------------------------------------------------------
    %-----------------------PROPAGATION/PREDICTION-----------------------------
    
    P_k_kminusone = F_d*P_kminusone_kminusone*F_d' + Q_d;
   %bibliography Maybeck cap 4.9

    
    %%-------------------------UPDATE/MEASUREMENT------------------------------
    
    % It is important to notice that we just get measurements in the velocity.
    % Then, we only write the H matrix for velocity measurements
    
    
    H = ConstructHv(L_hat,C_q_i_c_hat,C_q_w_i_hat,Cross_v_w_i_hat,Cross_omega_x,Cross_p_i_c_hat,p_i_c_hat,v_w_i_hat);
    
    
    %It is important noting that the resulting value from the projected
    %error state vector is an expecter measurement error. This means we
    %cannot take the difference ffrom this to the actual measurement, but
    %to the difference between the actual measurement and the expected
    %measurement value. 
    
    z_exp = L_hat *C_q_i_c_hat*(C_q_w_i_hat*v_w_i_hat' + Cross_omega_x*p_i_c_hat'); 
%     z_exp =0;
    %Compute the expected mean measurement.
    
    
    % Compute the residual in velocity
    r = vel_IMUOF(i,:)' - z_exp; 
    
    % Covariance Matrix of the Optical Flow measurements
    R; % defined in initialization of Kalman Filter
    
    %Compute the innovation
    S = H * P_k_kminusone * H' + R;
    
    % Compute Kalman Gain
    K  = P_k_kminusone*H'/S;
    
    % Compute correction
    X_correction = K*r;
    
    % Update Covariance Matrix
    P_k_k =  (eye(22) - K*H) * P_k_kminusone * (eye(22) - K*H)' + K*R*K';
    
    % Update Error State Vector
    X_err(i,:) = X_correction'; % again, the state is expressed horizontal
    
    %FIXME: this is a debug: Ground truth
    acc_gt(i,1:3) = IMUacc_gt(index_end,:);
    vel_gt(i,1:3) = vel(index_end+1,:);

    X_gt(i,:) = [ViconPose(index_end+3,1:3) vel_gt(i,1:3) ViconPose(index_end+3,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1];
        
    rotvel_gt = IMUrotvel(index_end,:);
    
    
    % Update State Vector
    X(i,1:6) = X_exp(i,1:6) + X_err(i,1:6);
    X(i,7:10) = quatmult(quat_normalize([0.5*X_err(i,7:9) 1]), X_exp(i,7:10));
    X(i,11:20) = X_exp(i,11:20) + X_err(i,10:19);
    X(i,21:24) = quatmult(quat_normalize([0.5*X_err(i,20:22) 1]), X_exp(i,21:24));
    
    
fprintf('Result of the Kalman Filter: \n');
fprintf('%d:X(i-1)_vel \t X_exp_vel \t Measurement \t Meas_exp \t X_corr \t X_post_vel \t X_gt_vel \n',i);

for ii=1:3
    fprintf('[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
             X(i-1,3+ii) ,X_exp(i,3+ii),vel_IMUOF(i,ii), z_exp(ii) ,  X_err(i,3+ii)  ,    X(i,3+ii) ,  X_gt(i,3+ii) );    
end

fprintf('%d:rot(i-1) \t X_exp_rot \t   -----  \t  ----- \t X_corr \t X_post_rot \t X_gt_rot \n',i);

for ii=1:3
    fprintf('[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
             X(i-1,6+ii) ,X_exp(i,6+ii),   0  , 0      ,  X_err(i,6+ii)  ,    X(i,6+ii) ,  X_gt(i,6+ii) );    
end
 fprintf('[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[ ---- ]\t[%6.5g]\t[%6.5g]\n',...
             X(i-1,10) ,X_exp(i, 10),  0  , 0,  X(i, 10) ,  X_gt(i, 10) );   

    % TODO: use IMU acceleration as a orientation measurement
   1;
   
end



%-------------------------Final Tests----------------------------------------
%% 3D trajectory

hPlot = figure('Name', sprintf('Position of Optical Flow'));
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
hold on
plot3(X(:,1),X(:,2),X(:,3),'r');
legend('Ground truth','EKF');
title('Optical Flow: Position');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis equal; grid on






%% Roll Pitch Yaw 
hPlot = figure('Name', sprintf('Orientation'));

[roll_EKF,pitch_EKF,yaw_EKF]=QuatToEuler(X(:,7:10));
[roll,pitch,yaw] = QuatToEuler(ViconPose(3:end,4:7));
roll = unwrap(roll);
pitch = unwrap(pitch);
yaw = unwrap(yaw);
% 
% roll_OF  =(POSE_OF(:,4));
% pitch_OF =(POSE_OF(:,5));
% yaw_OF  =(POSE_OF(:,6));
% 
% roll_IMUOF=(POSE_IMUOF(:,4));
% pitch_IMUOF=(POSE_IMUOF(:,5));
% yaw_IMUOF  =(POSE_IMUOF(:,6));
% 
% roll= (ViconPose(3:end,4));
% pitch =ViconPose(3:end,5);
% yaw = (ViconPose(3:end,4:6));



subplot(3,1,1);
title('Orientation')
plot(OpticalFlowTime(OpticalFlowstart:end),unwrap(roll_EKF(OpticalFlowstart:end)));
hold on;
plot(IMUTime,roll,'k')
xlabel('seconds'); ylabel('roll');
legend('EKF','Ground truth');
grid on

subplot(3,1,2);

plot(OpticalFlowTime(OpticalFlowstart:end),pitch_EKF(OpticalFlowstart:end));
hold on;
plot(IMUTime,pitch,'k')
xlabel('seconds'); ylabel('pitch');
legend('EKF','Ground truth');
grid on

subplot(3,1,3);

plot(OpticalFlowTime(OpticalFlowstart:end),yaw_EKF(OpticalFlowstart:end));
hold on;
plot(IMUTime,unwrap(yaw),'k')
xlabel('seconds'); ylabel('yaw');
legend('EKF','Ground truth');
grid on




%% Unrotated Linear Velocities

vel_truth =  diff(ViconPose(2:end,1:3))./timedifference(2:end,:);

hPlot = figure('Name', sprintf('Velocity of Optical Flow in world intertial frame'));
subplot(3,1,1);
title('Velocity of Optical Flow in world intertial frame');
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,1));grid on
hold on;
plot(IMUTime, vel_truth(:,1),'k');grid on
xlabel('seconds'); ylabel('velocity in x [mm/s]');
legend('EKF','Ground truth');

subplot(3,1,2);
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,2));grid on
hold on;
plot(IMUTime, vel_truth(:,2),'k');grid on
xlabel('seconds'); ylabel('velocity in y [mm/s]');
legend('EKF','Ground truth');
subplot(3,1,3);
plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,3));grid on
hold on;
plot(IMUTime, vel_truth(:,3),'k');grid on
xlabel('seconds'); ylabel('velocity in z [mm/s]');
legend('EKF','Ground truth');


%% Plot scale 

hPlot = figure('Name', sprintf('Scale Evolution'));
plot(OpticalFlowTime,X(:,17))


%% Rotation velocity

% hPlot = figure('Name', sprintf('Rotation Velocity'));
% 
% title('Rotation Velocity')
% subplot(3,1,1);
% plot(OpticalFlowTime(OpticalFlowstart:end),X(OpticalFlowstart:end,1)./OFTimedifference);grid on
% hold on;
% plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,1)./OFTimedifference,'r');grid on
% plot(IMUTime,IMUrotvel(:,1),'k');grid on
% xlabel('seconds'); ylabel('\omega_x [rad/s]');
% legend('5DoF','IMU and 3DoF','Ground truth');
% subplot(3,1,2);
% plot(OpticalFlowTime(OpticalFlowstart:end),omega5DoF(OpticalFlowstart:end,2)./OFTimedifference);grid on
% hold on;
% plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,2)./OFTimedifference,'r');grid on
% plot(IMUTime,IMUrotvel(:,2),'k');grid on
% xlabel('seconds'); ylabel('\omega_y [rad/s]');
% legend('5DoF','IMU and 3DoF','Ground truth');
% subplot(3,1,3);
% plot(OpticalFlowTime(OpticalFlowstart:end),omega5DoF(OpticalFlowstart:end,3)./OFTimedifference);grid on
% hold on;
% plot(OpticalFlowTime(OpticalFlowstart:end),omega_IMU(OpticalFlowstart:end,3)./OFTimedifference,'r');grid on
% plot(IMUTime,IMUrotvel(:,3),'k'); grid on
% xlabel('seconds'); ylabel('\omega_z [rad/s]');
% legend('5DoF','IMU and 3DoF','Ground truth');
% grid on
% 
% subplot(3,1,1);
% title('Rotation Velocity of Optical Flow')

