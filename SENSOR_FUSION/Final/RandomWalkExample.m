%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical Flow
% Cristian Duguet

clc
close all
clear all

% Gravitational Constant
G =  9.80665;

TRUE =1;
FALSE =0;

% Experiment name
experiment_date = '20120913_1545';
% experiment_date = '20120704_1913';
% experiment_date = '20120903_stones';

%% General Parameters:
%metrics in mm
% sigmaIMUacc = [0.6 0.6 0.9]*G*1e-3 * sqrt(IMUfreq); % [ in (mm/s2) / sqrt (Hz)  ]  * sqrt*(Hz)
% sigmaIMUrotvel = deg2rad([0.05 0.05 0.05]) * sqrt(IMUfreq); % as the noise specified in ADXRS610 [  (rad/s)/sqrt(Hz) ]  * sqrt(Hz)


% Focal Length (mm):
f = 6; %mm

% IMU Update Rate:
%IMU_rate = 50; %Hz

% Scaling constant to be changed
manual_scale = 36;

% Number of image instants which do not have Optical Flow
MAX_NOFLOW = 10;

% Pixel size
pixel_size=3.75e-3;


%% Read files

%If the mat files are not created from the log files, then create them
try
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_IMU.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
catch ME
    ReadFlowArrayFile(strcat(experiment_date,'_Flow.txt'));
    ReadIMU(strcat(experiment_date,'_IMU.txt'));
    ReadPoseVicon(strcat(experiment_date,'_PoseVicon.txt'));
    
    load(strcat(experiment_date,'_Flow.txt.mat'));
    load(strcat(experiment_date,'_IMU.txt.mat'));
    load(strcat(experiment_date,'_PoseVicon.txt.mat'));
end

%Analyze the frequency of messages
delta = diff(ViconTime,1);
Viconfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

delta = diff(OpticalFlowTime);
OpticalFlowfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));
% FIXME: sample rate of the camera does not seem to be 15FPS

delta = diff(IMUTime);
IMUfreq = mean(1./(double(delta(:,1)) + 1e-9*double(delta(:,2))));

fprintf('Mean Optical Flow rate: %g Hz\n',OpticalFlowfreq);
fprintf('Mean Vicon rate: %g Hz\n',Viconfreq);
fprintf('Mean IMU rate: %g Hz\n',IMUfreq);


% Coverting ground truth into mm
ViconPose(:,1:3) = ViconPose(:,1:3)*1000;

% Converting ground truth to initialize with orientation 0 0 0 1

ViconPose(:,4:7) = quatmult(ViconPose(:,4:7),repmat(invquat(ViconPose(1,4:7)),length(ViconTime),1));


% IMPORTANT: Converting angular velocities to rad/sec 
IMUrotvel = deg2rad(IMUrotvel);

% sigmaIMUacc = AngularVelocityCovariance;
% sigmaIMUrotvel = AccelerationCovariance; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Preparation

%-------------------------- Time adjustments ------------------------------

% We check and begin the experiment when the IMU and the IMU are delivering
% data. For that we will just 'cut off' the time when just one of the
% devices is broadcasting

ViconTime = double((ViconTime(:,1))) + double((ViconTime(:,2)))/1e9 ;
ViconSysTime = double((ViconSysTime(:,1))) + double((ViconSysTime(:,2)))/1e9 ;

OFTime = double((OpticalFlowTime(:,1))) + double((OpticalFlowTime(:,2)))/1e9;
OFSysTime = double((OpticalFlowSysTime(:,1))) + double((OpticalFlowSysTime(:,2)))/1e9;

IMUTime = double((IMUTime(:,1))) + double((IMUTime(:,2)))/1e9 ;
IMUSysTime = double((IMUSysTime(:,1))) + double((IMUSysTime(:,2)))/1e9 ;


%FIXME : JUST FOR THIS EXPERIMENT 
%The clock with which IMUTime was saved was much more delayed. We have to
%correct this clock bias assuming that the delay of the IMU signals were
%about the same 
OF_delay_avg = mean(OFSysTime-OFTime);
IMU_delay_avg = mean(IMUSysTime-IMUTime);
Vicon_delay_avg = mean(ViconSysTime-ViconTime);

IMUTime = IMUTime + IMU_delay_avg - OF_delay_avg;
ViconTime = ViconTime + Vicon_delay_avg - OF_delay_avg;


% the time in which both are transmiting
time_zero =  min([OFTime(1) IMUTime(1) ViconTime(1)]);
% 
OFTime = OFTime-time_zero;
OFSysTime = OFSysTime- time_zero;
IMUTime = IMUTime-time_zero;
IMUSysTime - IMUSysTime - time_zero;
ViconTime = ViconTime-time_zero;
ViconSysTime = ViconSysTime-time_zero;



%-----------------Simulated IMU Construction ------------------------------


[SIMUacc, SIMUrotvel,SIMUTime] = Construct_SIMU(ViconTime,ViconPose,AccelerationCovariance,AngularVelocityCovariance);
SIMUSysTime = ViconSysTime(3:end);




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

% For Optical Flow
noflow_counter =0;
max_inliers = zeros(size(matched_feat));

% % Initial Pose (Initial state Asumption) in the form (X,Y,Z,rx,ry,rz,w)
% POSE_OF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison
% POSE_IMUOF = repmat([ViconPose(1,1:3) 0 0 0 1],length(OpticalFlowTime),1);  %FIXME: We are taking the ground truth initial position, for comparison

%Sensors Information initialization
vel_IMUOF = zeros(length(OpticalFlowTime),3);
omega_IMU = zeros(length(OpticalFlowTime),3);

X_gt = repmat([ViconPose(1,1:3) 0 0 0 ViconPose(1,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1],length(OpticalFlowTime),1);
vel_gt = repmat( [0 0 0],length(OpticalFlowTime),1);
acc_gt = repmat( [0 0 0],length(OpticalFlowTime),1);

omega_m = zeros(length(OpticalFlowTime),3);
a_m = zeros(length(OpticalFlowTime),3);


%% ------------------KALMAN FILTER INITIALIZATION -------------------------

% we  express the first state as X0, the prior as X_exp, the
% posterior state as X

%------------------------------Noises--------------------------------------

EKF = Indirect_ExtendedKalmanFilter;

%Acceleration noise variance mm/s2
EKF.var_na = [3 3 3];
% Acceleration Bias noise variance mm/s2
sigma_nba = [0.01 0.01 0.01];
EKF.var_nba = sigma_nba.^2;
% sigma_nba = [0 0 0];
% Rotational velocity noise variance
EKF.var_nomega = [1 1 1];
%Rotational velocity bias noise variance
sigma_nbomega = [0.1 0.1 0.1];   
EKF.var_nbomega = sigma_nbomega.^2;

%scale noise
EKF.var_nL = 0.05^2; % relative to unit

% Covariance Matrix of the Optical Flow measurements
EKF.R = diag([40 40 90].^2); %mm2
%FIXME: better representation

% Noises needed by Stephan weiss q_d calculationa lgorithm
EKF.var_nqvw = [0 0 0];
EKF.var_nqci = [0 0 0];
EKF.var_npic = [0 0 0];

EKF.Image_width = Image_width ;
EKF.Image_height = Image_height ;

%-----------------------State Initialization-------------------------------

%We assume the initial position is at one meter height
p_w_i   = ViconPose(1,1:3);   %FIXME: We are taking the ground truth initial position, for comparison
v_w_i   = [0 0 0]; % Assume zero velocity  %FIXME
q_w_i   = [0 0 0 1];  % Assume no rotation at the beginning
b_omega = [0 0 0];
b_a     = [0 0 0];
L       = 0.005;
p_i_c   = [0 -4 -4];   %Calibration state: Distance from IMU to camera
q_i_c   = [0 0 0 1];     %calibration state: Orientation of camera w.r.t IMU frame

% For the initial state, the state vector wil be
EKF.X = [ p_w_i v_w_i q_w_i b_omega b_a L p_i_c q_i_c] ;  % 24 element vector
%Expected state vector
% For the initial state, the error state covariance matrix will be
% trueerror = zeros(1,22);
% trueerror(1:6) = X0(1:6)-X_gt(1:6);
% aux = quatmult(X0(7:10),invquat(X_gt(6:10)));
% trueerror(7:9) = aux(1:3);
% trueerror(10:19)=X0(11:20)-X_gt(11:20);
% aux = quatmult(X0(21:24),invquat(X_gt(21:24)));
% trueerror(20:22) = aux(1:3);
% P_kminusone_kminusone = trueerror'*trueerror +.5*ones(22);
% P_kminusone_kminusone = 0.1*eye(22);

EKF.P_k_k = [0.016580786012789, 0.012199934386656, -0.001458808893504, 0.021111179657363, 0.007427567799788, 0.000037801439852, 0.001171469788518, -0.001169015812942, 0.000103349776558, -0.000003813309102, 0.000015542937454, -0.000004252270155, -0.000344432741256, -0.000188322508425, -0.000003798930056, 0.002878474013131, 0.000479648737527, 0.000160244196007, 0.000012449379372, -0.000025211583296, -0.000029240408089, -0.000001069329869, -0.001271299967766, -0.000133670678392, -0.003059838896447
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
EKF.P_k_k = EKF.P_k_k(1:22,1:22);

%FIXME




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin of the iterative algorithm


% BIG ASSUMPTION: WE ASSUME HERE THAT PROCESSING TIME IS 0 . SO THIS MEANS
% WE ARE JUST CONSIDERING THE SYSTEM DELAYS OUTSIDE OF THE FILTER.
% (EKF TIME, 8POINT TIME CONSIDERED ZERO)

% SECOND ASSUMPTION: MESSAGES OF EACH TOPIC ARE ARRIVING IN CONSECUTIVE
% FORM , I.E. MEASUREMENT K ALWAYS ARRIVE BEFORE MEASUREMENT K+1


% 
TimeSteps_OF = length(OFTime);
TimeSteps_SIMU = length(SIMUTime);

k_OF = 2;
k_SIMU = 2;

%%-------------------------------------------------------------------------
% Save the Filter result in the following vectors 
% FFilter time: 
t_EKF = nan(TimeSteps_OF + TimeSteps_SIMU - 2,1);
k_EKF = 1;
X =  repmat([ p_w_i v_w_i q_w_i b_omega b_a L p_i_c q_i_c],length(t_EKF),1);
%%-------------------------------------------------------------------------

thres_delay = 100e-3; % 100 ms threshold tolerance 

tic % Begin clock count 
% for k_SIMU =2:TimeSteps_SIMU
%      EKF = EKF.Propagate(SIMUrotvel(k_SIMU-1:k_SIMU,:), SIMUacc(k_SIMU-1:k_SIMU,:), SIMUTime(k_SIMU) - SIMUTime(k_SIMU-1));
%         X(k_EKF,:) = EKF.X_exp;
%         EKF.X = EKF.X_exp;
%         t_EKF(k_EKF) = SIMUTime(k_SIMU);
%         k_EKF = k_EKF+1;
%         t_real = SIMUSysTime(k_SIMU);
%         EKF.propagated = 1;
% end
% 
% 

while(k_OF<=TimeSteps_OF && k_SIMU <=TimeSteps_SIMU)
    if OFSysTime(k_OF) >SIMUSysTime(k_SIMU)
        % if the message corresponds to a past measurement, and it is
        % already too old, discard it
        if SIMUTime(k_SIMU) < t_EKF
            k_SIMU = k_SIMU+1;
            continue
        end
        
        % Wait the treshold delay time and if something "earlier" comes,
        % then update first instead of propagate
        if OFSysTime(k_OF) - SIMUSysTime(k_SIMU) > thres_delay && OFTime(k_OF) < SIMUTime(k_SIMU)
%             fprintf('updating...\n');
%             EKF = EKF.Update(OpticalFlow_x(k_OF,1:matched_feat(k_OF)), OpticalFlow_y(k_OF,1:matched_feat(k_OF)), ...
%                 OpticalFlow_u(k_OF,1:matched_feat(k_OF)), OpticalFlow_v(k_OF,1:matched_feat(k_OF)))    ;
%             X(k_EKF,:) = EKF.X;
%             t_EKF(k_EKF) = OFTime(k_OF);
%             k_EKF = k_EKF+1;
%             t_real = max(t_real,OFSysTime(k_SIMU));
            k_OF = k_OF +1;
            EKF.propagated = 0;
            continue
        end
        
        fprintf('propagating...\n');
        EKF = EKF.Propagate(SIMUrotvel(k_SIMU-1:k_SIMU,:), SIMUacc(k_SIMU-1:k_SIMU,:), SIMUTime(k_SIMU) - SIMUTime(k_SIMU-1));
        X(k_EKF,:) = EKF.X_exp;
        EKF.X = EKF.X_exp;
        t_EKF(k_EKF) = SIMUTime(k_SIMU);
        k_EKF = k_EKF+1;
        t_real = SIMUSysTime(k_SIMU);
        k_SIMU = k_SIMU+1;
        EKF.propagated = 1;
        
    elseif SIMUSysTime(k_SIMU) >= OFSysTime(k_OF)  && EKF.propagated == 1
        
        
        % if the message corresponds to a past measurement, and it is
        % already too old, discard it
        if OFTime(k_OF) < t_EKF
            k_OF = k_OF+1;
            continue
        end
        
        % Wait the treshold delay time and if something "earlier" comes,
        % then propagate first instead of update
        if SIMUSysTime(k_SIMU) - OFSysTime(k_OF) > thres_delay && OFTime(k_OF) > SIMUTime(k_SIMU)
            
            fprintf('propagating...\n');
            EKF = EKF.Propagate(SIMUrotvel(k_SIMU-1:k_SIMU,:), SIMUacc(k_SIMU-1:k_SIMU,:), SIMUTime(k_SIMU) - SIMUTime(k_SIMU-1));
            X(k_EKF,:) = EKF.Xexp;
            EKF.X = EKF.X_exp;
            t_EKF(k_EKF) = SIMUTime(k_SIMU);
            k_EKF = k_EKF+1;
            t_real = max(t_real,SIMUSysTime(k_SIMU));
            k_SIMU = k_SIMU+1;
            EKF.propagated = 1;
            continue
        end
        fprintf('updating...\n');
%         EKF = EKF.Update(OpticalFlow_x(k_OF,1:matched_feat(k_OF)), OpticalFlow_y(k_OF,1:matched_feat(k_OF)), ...
%                 OpticalFlow_u(k_OF,1:matched_feat(k_OF)), OpticalFlow_v(k_OF,1:matched_feat(k_OF)))    ;
%         X(k_EKF,:) = EKF.X;
%         t_EKF(k_EKF) = OFTime(k_OF);
%         k_EKF = k_EKF+1;
%         t_real = max(t_real,OFSysTime(k_OF));
        k_OF = k_OF +1;
        EKF.propagated = 0;

    else  %in case we haven't propagated yer we skip the update
        k_OF = k_OF+1;
    end
    
    
end

% TODO: use IMU acceleration as a orientation measurement
fprintf('Finished Experiment. \n');
fprintf('Experiment real time elapsed: %d \t Simulation Time elapsed: %d \n',t_real-SIMUSysTime(2),toc);
fprintf('Result of the Kalman Filter: \n');

% fprintf('%d:Xprev_vel \t X_exp_vel \t Measurement \t Meas_exp \t X_corr \t X_post_vel \t X_gt_vel \n',i);
% 
% for ii=1:3
%     fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
%         X(i-1,3+ii) ,X_exp(i,3+ii),vel_IMUOF(jj,ii), z_exp(ii) ,  X_err(i,3+ii)  ,    X(i,3+ii) ,  X_gt(i,3+ii) );
% end
% 
% fprintf('%d:rot(i-1) \t X_exp_rot \t   -----  \t  ----- \t X_corr \t X_post_rot \t X_gt_rot \n',i);
% 
% for ii=1:3
%     fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
%         X(i-1,6+ii) ,X_exp(i,6+ii),   0  , 0      ,  X_err(i,6+ii)  ,    X(i,6+ii) ,  X_gt(i,6+ii) );
% end
% fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[ ---- ]\t[%6.5g]\t[%6.5g]\n',...
%     X(i-1,10) ,X_exp(i, 10),  0  , 0,  X(i, 10) ,  X_gt(i, 10) );
% 
% 

fprintf('Done!')

%         
%         %------------------------------------------------------------------
%         
%         
%         
%         
%         
%         
%         %FIXME: this is a debug: Ground truth
%         acc_gt(i,1:3) = SIMUacc_no_g(i,:);
%         vel_gt(i,1:3) = vel(i+1,:);
%         
%         X_gt(i,:) = [ViconPose(i+3,1:3) vel_gt(i,1:3) ViconPose(i+3,4:7) 0 0 0 0 0 0 1 -40 -40 -40 0 0 0 1];
%         
%         rotvel_gt = IMUrotvel(i,:);
%         
%         
%        
%         
%         % TODO: use IMU acceleration as a orientation measurement
%         fprintf('Result of the Kalman Filter: \n');
%         fprintf('%d:Xprev_vel \t X_exp_vel \t Measurement \t Meas_exp \t X_corr \t X_post_vel \t X_gt_vel \n',i);
%         
%         for ii=1:3
%             fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
%                 X(i-1,3+ii) ,X_exp(i,3+ii),vel_IMUOF(jj,ii), z_exp(ii) ,  X_err(i,3+ii)  ,    X(i,3+ii) ,  X_gt(i,3+ii) );
%         end
%         
%         fprintf('%d:rot(i-1) \t X_exp_rot \t   -----  \t  ----- \t X_corr \t X_post_rot \t X_gt_rot \n',i);
%         
%         for ii=1:3
%             fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\n',...
%                 X(i-1,6+ii) ,X_exp(i,6+ii),   0  , 0      ,  X_err(i,6+ii)  ,    X(i,6+ii) ,  X_gt(i,6+ii) );
%         end
%         fprintf(' [%6.5g]\t[%6.5g]\t[%6.5g]\t[%6.5g]\t[ ---- ]\t[%6.5g]\t[%6.5g]\n',...
%             X(i-1,10) ,X_exp(i, 10),  0  , 0,  X(i, 10) ,  X_gt(i, 10) );
%         
%         continue
% end
    
%     X(i,:) = X_exp(i,:);
%     P_kminusone_kminusone = P_k_kminusone;
 
    % TODO: use IMU acceleration as a orientation measurement
   
    


%-------------------------Final Tests----------------------------------------
%% 3D trajectory

hPlot = figure('Name', sprintf('Position of EKF'));
plot3(ViconPose(:,1),ViconPose(:,2),ViconPose(:,3));
hold on
plot3(X(:,1),X(:,2),X(:,3),'r');
legend('Ground truth','EKF');
title('Optical Flow: Position');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis equal; grid on

% 3D Trajetory of Optical Flow 
% hPlot = figure('Name', sprintf('Position of Optical Flow'));
% 
% plot3(1/L_hat*cumsum(vel_IMUOF(:,1)),1/L_hat*cumsum(vel_IMUOF(:,2)),1/L_hat*cumsum(vel_IMUOF(:,3)))
% xlabel('x');ylabel('y');zlabel('z');

%% Roll Pitch Yaw
hPlot = figure('Name', sprintf('Orientation'));

[roll_EKF,pitch_EKF,yaw_EKF]=QuatToEuler(X(1:k_EKF,7:10));
[roll,pitch,yaw] = QuatToEuler(ViconPose(:,4:7));
roll = unwrap(roll);
pitch = unwrap(pitch);
yaw = unwrap(yaw);

[roll_SIMU,pitch_SIMU,yaw_SIMU]=QuatToEuler(IMUrot);



% roll_OF  =(POSE_OF(:,4));
% pitch_OF =(POSE_OF(:,5));
% yaw_OF  =(POSE_OF(:,6));
%
% roll_SIMUOF=(POSE_SIMUOF(:,4));
% pitch_SIMUOF=(POSE_SIMUOF(:,5));
% yaw_SIMUOF  =(POSE_SIMUOF(:,6));

% roll= (ViconPose(3:end,4));
% pitch =ViconPose(3:end,5);
% yaw = (ViconPose(3:end,4:6));


subplot(3,1,1);
title('Orientation')
plot(t_EKF(1:k_EKF),unwrap(roll_EKF));
hold on;
plot(ViconTime,roll,'k')
plot(IMUTime,roll_SIMU,'r');
xlabel('seconds'); ylabel('roll');
legend('EKF','Ground truth','SIMUfilter');
grid on

subplot(3,1,2);

plot(t_EKF(1:k_EKF),pitch_EKF);
hold on;
plot(ViconTime,pitch,'k')
plot(IMUTime,yaw_SIMU,'r');
xlabel('seconds'); ylabel('pitch');
legend('EKF','Ground truth','SIMUfilter');
grid on

subplot(3,1,3);

plot(t_EKF(1:k_EKF),yaw_EKF);
hold on;
plot(ViconTime,unwrap(yaw),'k')
plot(IMUTime,yaw_SIMU,'r');
xlabel('seconds'); ylabel('yaw');
legend('EKF','Ground truth','SIMUfilter');
grid on




%% Unrotated Linear Velocities

vel_truth =  diff(ViconPose(:,1:3))./repmat(diff(ViconTime),1,3);

hPlot = figure('Name', sprintf('Velocity of Optical Flow in world intertial frame'));
subplot(3,1,1);
title('Velocity of Optical Flow in world intertial frame');
plot(t_EKF,X(:,4));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,1),'k');grid on
% plot(SIMUTime,SIMU(:,),'r');
xlabel('seconds'); ylabel('velocity in x [mm/s]');
legend('EKF','Ground truth');

subplot(3,1,2);
plot(t_EKF,X(:,5));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,2),'k');grid on
% plot(OpticalFlowTime,1/L_hat*vel_SIMUOF(:,1),'-g');
xlabel('seconds'); ylabel('velocity in y [mm/s]');
legend('EKF','Ground truth');
subplot(3,1,3);
plot(t_EKF,X(:,6));grid on
hold on;
plot(ViconTime(2:end), vel_truth(:,3),'k');grid on
% plot(OpticalFlowTime,1/L_hat*vel_SIMUOF(:,1),'-g');
xlabel('seconds'); ylabel('velocity in z [mm/s]');
legend('EKF','Ground truth');


%% Plot scale

hPlot = figure('Name', sprintf('Scale Evolution'));
plot(IMUTime,X(:,17));


% Plot rotation bias

hPlot = figure('Name', sprintf('Rotation velocity bias'));
title('Rotation Velocity Bias');
plot(IMUTime,X(:,11:13))
legend('x','y','z');

%plot acceleration bias

hPlot = figure('Name', sprintf('acceleration bias'));
title('Acceleration Bias');
plot(IMUTime,X(:,14:16))
legend('x','y','z');


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

