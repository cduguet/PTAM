% Cristian Duguet
% Particle Filter

classdef ParticleFilter
    properties
        
        % States have the structure X = [ p_w_i v_w_i q_w_i b_omega b_a]
        
        %==================PARTICLE PARAMETERS=============================
        
        % Number of Particles
        N_particles;
        
        %effective number of particles
        N_eff;
        
        % Particles array
        X;
        
        % Weights
        weights;
        
        %=====================NOISE AND STANDARD DEVIATIONS================
        
        % Acceleration noise std dev
        sigma_na;
        % Acceleration bias noise std dev
        sigma_nba;
        % Rotation velocity std dev
        sigma_nomega;
        % Rotation velocity std dev 
        sigma_nbomega;
        %Acceleration bias noise variance
        %         var_nba;
        %angular velocity variance
        %         var_nomega;
        %angular velocity bias variance
        %         var_nbomega;
        
        % Optical Flow covariance matrix
        sigma_OF;
        
        %-------------------------- state parameters-----------------------
        
        % Scale of  the Optical Flow
        L =1;
        
        % position of camera w.r.t. to imu frame
        p_i_c;
        % rotation of camera w.r.t. to imu frame
        q_i_c;
        
        %Extra state from last update
        q_w_i_lastOF ;        
        
        z_OF = [0 0 0]';
        
        %------------------------fIR low pass filter-----------------------
        
        FIR_buffsize =5;
        FIR_coeff = [.2 .2 .2 .2 .2];        
        FIR_buff;
                
        %---------------------- counters and flags-------------------------
        % Propagated flag
        propagated = 0;
        
        % number of resamplings
        resample_counter =0;
        
        %tolerance to optical flow
        noflow_counter =0;
        
        %----------------------- Optical Flow Parameters--------------------
        
        % Size of the pixel
        pixel_size = 3.75e-6;
        
        % Focal length
        f = 6e-3;
        
        %Image size
        Image_height;
        Image_width;
        
        %---------------------------- Ransac Parameters--------------------
        
        max_Iter = 100;
        sample_num = 50;
        RANSAC_Linear_treshold = 1;
        %---------------------------Debug variables------------------------
        z_SOF; 
        
        propagatetime=0;
        opticalflowsimtime=0;
        outlierstime=0;
        updatetime=0;
        resampletime=0;
        %------------------------------------------------------------------
        
    end
    
    properties (Dependent = true)
        %Weighted mean of Posterior
        X_mean;
    end
    
    methods
        %% PROPAGATE
        function PF = Propagate(PF, gyro, acc, delta_t)
            %% ================= IMU Strapdown Algorithm ==========================
            % The initial velocity and position information is in X( i-1,:)
            % INPUT THE LAST TWO MEASUREMENTS, SINCE IT IS INTEGRATING IN
            % 1ST ORDER
            tic
            for i = 1:PF.N_particles
                PF.X(i,:) = IMU_PFStrapdown(PF.X(i,:)...
                    ,gyro+ [0 0 0 ; randn(1,3) * diag(PF.sigma_nomega)]...
                    ,acc + [0 0 0 ; randn(1,3) * diag(PF.sigma_na)]...
                    ,delta_t);
            end
            
            PF.X(:,11:13) = PF.X(:,11:13) + randn(PF.N_particles,3) * diag(PF.sigma_nbomega);
            PF.X(:,14:16) = PF.X(:,14:16) + randn(PF.N_particles,3) * diag(PF.sigma_nba);
            timeaux=toc;
            PF.propagatetime =PF.propagatetime +timeaux;
        end
        
        
        %% UPDATE
        
        function PF = Update(PF, x, y, u, v, delta_t)
%             tic
            % --------------------- Expected State Values -----------------
            X_prior = PF.X_mean;
            q_w_i       = X_prior(7:10);  % Expectation fof the quaternion from world to IMU
            b_omega     = X_prior(11:13); % Expectation of the bias of the IMU gyro measurements
            
            % Mean rotation being made between updates
            Rotation = quatmult(q_w_i,invquat(PF.q_w_i_lastOF));
            % convert this rotation to a mean angular velocity
            omega_mean  = 2* Rotation(1:3)/ delta_t;
            
            Rotation = QuatToRotMat(Rotation);
            % update last rotation to it's new value
            PF.q_w_i_lastOF = q_w_i;
            
            %calulate cross matrix from mean rotation velocity
            Cross_omega_x = VectorToCrossMat(omega_mean - b_omega);
            
            
            %Make rotation matrices from quaternions
            C_q_w_i = QuatToRotMat(q_w_i);
            C_q_i_c = QuatToRotMat(PF.q_i_c);
            
            
            
            %% =====================OPTICAL FLOW===========================
            
            ama=PF.pixel_size; % FIXME : get rid of this constant
            
            % Set the pixel position referenced to the center
            x = x'-PF.Image_width /2+0.5;
            y = y'-PF.Image_height/2+0.5;
            x = x *(ama/PF.f); y = -y *(ama/PF.f);
            z = ones(size(x));
            
            u = u'*(ama/PF.f);
            v = -v'*(ama/PF.f);
               
            % In case of no flow detected
            if length(u) < 5
                %                if     noflow_counter > MAX_NOFLOW
                %                    error('Optical Flow information too poor. Need more features.\n');
                %                    return
                %                end
                %                noflow_counter  = noflow_counter+1;
                warning('No flow detected in this image frame.\n');
                return
            end
            %             toc1 = toc;
            
            %% ---------------------OUTLIERS REJECTION---------------------
            % ----------------------Homography RANSAC----------------------
            % First method: Calculating the Homogrpahy Matrix from a random subset
            % of sample_num points of all the matched features in the optical flow.
            % After max_Iter iterations, return the inliers of the most accepted
            % model (The model which can include the most inliers).
            
            % inliers_idx = RANSAC_Homography(x,y,u,v,max_Iter,sample_num,RANSAC_Homography_treshold);
            %----------------------Linear RANSAC---------------------------
            % Second method: Matching the Consensus of the random samples in a
            % linear model
            % tic
            
            % inliers_idx = RANSAC_Linear(u,v,PF.max_Iter,PF.sample_num,PF.RANSAC_Linear_treshold);
            
            %--------------------- Least Medaian of Squares Rejection------
            tic
            
            
            %---------------------Histogram Analysis-----------------------
            inliers_idx = HistographicRejection(u,v);
            
            
            %--------------------------------------------------------------
            
            if length(inliers_idx) < 5
                warning('Inliers information too poor (%d). Avoiding outlier rejection.\n',length(inliers_idx));
                w = zeros(size(u));
            else
                %                full_x=x;
                %                full_y=y;
                %                full_u=u;
                %                full_v=v;
                x = x(inliers_idx);
                y = y(inliers_idx);
                z = z(inliers_idx);
                u = u(inliers_idx);
                v = v(inliers_idx);
                w = zeros(size(u));
            end
            
            
            timeaux = toc;
            PF.outlierstime = PF.outlierstime + timeaux;
            
            %% ------------ SENSOR COLLIGATION ----------------------------
            tic
            % Unrotate Optical Flow for use in 2DoF Epipolar Algorithm
            % The rotation velocity of the IMU has been integrated over time. The
            % result of it is the rotation quaternion. As the Optical Flow here is
            % a distance vector rather than a velocity vector, we use the
            % integrated rotation of the IMU to unrotate the Optical Flow
            % Assuming smallangle:
            %             tic
            
            unrotation_OF = (eye(3)-Rotation) * [x y z]';
            
            %UNROTATE
            u = u + unrotation_OF(1,:)';
            v = v + unrotation_OF(2,:)';
            w = w + unrotation_OF(3,:)';
            
            %% --------------and EPIPOLAR RECONSTRUCTION ------------------
            
            %Calculate using the Epipolar Constraint over the Optical Flow
            
            %------------------ planar Optical Flow (w=0)------------------
%             aux = Continuous8point1DoF(x,y,z,u,v,w)/delta_t;
            
            %This option demonstrated to be the best: 
            
            aux = Continuous8point2DoF(x,y,z,u,v,w)/delta_t;
            %----------------FIR Filter------------------------------------
            
            NLdistort = 2*1/PF.L;
            
            aux = PF.FIR_buff(:,1) ...
                + sign(((aux - PF.FIR_buff(:,1))/2)).*(1 -exp(NLdistort*-abs(-aux + PF.FIR_buff(:,1))))/NLdistort;
            
            
            PF.FIR_buff = [aux PF.FIR_buff(:,1:4)];
            aux = (PF.FIR_buff * PF.FIR_coeff');
            
                                    
            
%             if norm(aux - PF.z_OF) > 2
%                 PF.z_OF = mean([PF.z_OF aux],2);
%             else
%                 PF.z_OF = aux;
%             end
            
            % ----------------Update---------------------------------------
            
            % Transform Optical Flow into metric measurements

            PF.z_OF = C_q_w_i'*C_q_i_c'*(PF.L* aux - C_q_i_c*Cross_omega_x *PF.p_i_c);
            PF.opticalflowsimtime = PF.opticalflowsimtime +toc;
            
            fprintf('OF meas: [%d %d %d]\n', PF.z_OF)
            
            tic
            %Rotate Covariance of the measurements
            Z_cov = (diag(C_q_w_i'*C_q_i_c'*PF.sigma_OF))^2;
            
            %Update the weights
            PF.weights = PF.weights.* mvnpdf(PF.X(:,4:6), repmat(PF.z_OF',PF.N_particles,1) , repmat(Z_cov,[1,1,PF.N_particles]));

            % Update  weights using x-y
%             PF.weights = PF.weights.* mvnpdf(PF.X(:,4:5), repmat(PF.z_OF(1:2)',PF.N_particles,1) , repmat(Z_cov(1:2,1:2),[1,1,PF.N_particles]));

            %%--debugging code---------------------------------------------
%             size(PF.z_SOF)
%             size(PF.X(3))
%             size(aux)
%             size(Cross_omega_x)
%             size(PF.p_i_c)
            

%             size(Z_cov)
            %Update using SOF
%             PF.z_SOF = C_q_w_i'*C_q_i_c'*(PF.z_SOF' - C_q_i_c*Cross_omega_x *PF.p_i_c);
%             PF.weights = PF.weights.* mvnpdf(PF.X(:,4:6), repmat(PF.z_SOF',PF.N_particles,1) , repmat(Z_cov,[1,1,PF.N_particles]));
%             if isnan(PF.z_SOF)
%                 error('OF measurement is NAN');
%             end
            
            %%--debugging code-end-----------------------------------------
            
            %% Draw 3D point cloud
%             hPlot = figure('Name', 'Particle Filter Point Cloud')
            
%             h1 = scatter3(PF.X(:,4),PF.X(:,5),PF.X(:,6),2,'b');
%             hold on;
%             h2 = scatter3(PF.z_OF(1),PF.z_OF(2),PF.z_OF(3),5,'r');
%             h3 = scatter3(PF.z_SOF(1),PF.z_SOF(2),PF.z_SOF(3),5,'*k');
%             hold off;
%             legend([h1 h2 h3],'Prediction','OF_meas','SOF_meas');
%             xlabel('x (m/s)');ylabel('y (m/s)');zlabel('z (m/s)');
            
%             
%             hold on;
%             for a = 10.^(1:4) % 'a' defines the isosurface limits
%                 p = patch(isosurface(X,Y,Z,V,max(max(max(V)))/a)); % isosurfaces at max(V)/a
%                 isonormals(X,Y,Z,V,p); % plot the surfaces
%                 set(p,'FaceColor','red','EdgeColor','none'); % set colors
%             end
%             alpha(.1); % set the transparency for the isosurfaces
%             daspect([1 1 1]); box on; axis tight;
            
            %%
            if (sum(PF.weights)) ==0
                error('weights sum is 0! \n');
            end
            
            %normalize weights
            PF.weights  = PF.weights/sum(PF.weights);
            
            %Calculate the effective Number of Particles
            PF.N_eff = 1/sum(PF.weights.^2);
            
            PF.updatetime = PF.updatetime +toc;

            if PF.N_eff < 0.5*PF.N_particles;
                tic
                fprintf('Resampling...\n');
                PF.Resample2;
                PF.resampletime = PF.resampletime + toc;
            end
            
            if isnan(sum(PF.weights))
                error('weights are NaN! \n');
            end
            
            
        end
        
        function Resample(PF)
            %cumulative function of weights
            cum = cumsum(PF.weights);
            %uniformely distributed intervals
            u = linspace(0,1-1/PF.N_particles,PF.N_particles) + unifrnd(0, 1/PF.N_particles);
            
            %initialize value asignation matrices
            X_P_new = zeros(size(PF.X));
            w_k_new = zeros(PF.N_particles,1);
            for h = 1:PF.N_particles
                index_found = find(cum>u(h),1);
                X_P_new(h,:) = PF.X(index_found,:);
                w_k_new(h) = PF.weights(index_found);
            end
            PF.X = X_P_new;
            PF.weights = w_k_new;
            PF.resample_counter =PF.resample_counter+1;
        end
        
        
        function Resample2(PF)
            
            index_max = find(PF.weights == max(PF.weights),1);
            PF.X = repmat(PF.X(index_max,:),length(PF.weights),1);
            PF.weights = repmat(1/PF.N_particles,length(PF.weights),1);
%             %cumulative function of weights
%             cum = cumsum(PF.weights);
%             %uniformely distributed intervals
%             u = linspace(0,1-1/PF.N_particles,PF.N_particles) + unifrnd(0, 1/PF.N_particles);
%             
%             %initialize value asignation matrices
%             X_P_new = zeros(size(PF.X));
%             w_k_new = zeros(PF.N_particles,1);
%             for h = 1:PF.N_particles
%                 index_found = find(cum>u(h),1);
%                 X_P_new(h,:) = PF.X(index_found,:);
%                 w_k_new(h) = PF.weights(index_found);
%             end
%             PF.X = X_P_new;
%             PF.weights = w_k_new;
%             PF.resample_counter =PF.resample_counter+1;
        end
        
        %------------------------------------------------------------------
        %%  OTHER
        %calculate mean of particles
        function X_mean = get.X_mean(PF)
            X_mean = zeros(1,16);
            %Mean of the linear components:
            X_mean(1:6) = sum(repmat(PF.weights,1,6).*PF.X(:,1:6),1);
            X_mean(11:16) = sum(repmat(PF.weights,1,6).*PF.X(:,11:16),1);
            
            %mean of Quaternions
            X_mean(7:10) = QuatAvg(PF.X(:,7:10),PF.weights);
        end
    end
end