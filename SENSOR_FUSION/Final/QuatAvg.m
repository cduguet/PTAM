% Cristian Duguet
% weighted Mean of Quaternions 

function  q_avg =  QuatAvg(Q,w)
%Make square mattrix qhich is the sum of the multiplication of the vertical
%quaternion q * q'. As we stored the quaternions in an horizontal fasion
%here, we calculate (including also the weights). \cite{markley2007}

M = Q'*diag(w)*Q;

% QUEST
%Quaternion Estimation Algorithm 
%To solve using Quest, use  matrix K = 4*M - sum(w) * eye(4)
% sum(w) is supposed to be 1
% reference to the problem in \cite{Hall2003} Chapter 4-18 

K = 4*M - sum(w) * eye(4);
sigma = K(4,4);
lambda_opt = sum(w);
S = K(1:3,1:3) + sigma * eye(3);
%Obtain the Rodriguez parameters
p = ((lambda_opt + sigma)*eye(3) - S) \ K(1:3,4);

%Construct the Optimal vector out of the Rodriguez parameters
q_avg = 1/sqrt(1+p'*p) *  [p;1];

if isnan(sum(q_avg))
    fprintf('Rodriguez parameters became singular. Calculating linear mean... \n');
    q_avg = sum(diag(w)*Q,1);
    
end
q_avg = quat_normalize(q_avg);
end 