
Main code of the Master's thesis. 
Cristian Duguet 

I implemented two codes for monocular and inertial sensor fusion to obtain velocity estimation. Real experiments data was saved with ROSbag and the day of the newest experiments is 2012-09-27. Rosbag files saved: 
imu data, 
distortion corrected image sequence, 
camera inforamtion, 
vicon position. 

The code in this folder takes the saved simulated data, and implements an Kalman and a Particle Fitler to estimate the state composed mainly of  pose, velocity, rotation velocity, and imu sensors' biases. 


for the EKF Filter, run the script Main_EKF.m in Matlab

for the Particle Filter, run the script Main_ParticleFilter.m in Matlab 