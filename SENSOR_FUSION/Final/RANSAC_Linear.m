% Calculate the inliers and outliers of a group of optical Flow points
% using an Homography model

function max_inliers_idx = RANSAC_Linear(u,v,maxIter,sample_num,RANSAC_treshold)

max_inliers =0;
matched_feat = length(u);
max_inliers_idx = [];

for j = 1:maxIter
            randompos = ceil(single(matched_feat)* rand(1,sample_num));
            random_u = u(randompos);
            random_v = v(randompos);
            
      %% Estimate the mean and standard deviations of the set
      mean_u = mean(random_u);
      mean_v = mean(random_v);
      sigma2_u = var(u);
      sigma2_v = var(v);
        
      %% Discern between out/ inliers 
      
      inliers_idx = find( (u-mean_u).^2 + (v-mean_v).^2  < (sigma2_u + sigma2_v));
      inliers_no = length(inliers_idx);
      if max_inliers < inliers_no
          max_inliers = inliers_no;
          max_mean_u = mean_u;
          max_mean_v = mean_v;
          max_inliers_idx = inliers_idx;
      end
      
%        proboutlier  = 1- max_inliers/double(matched_feat);
%        e = min(e,proboutlier);
%        loopno = log10(1-p)/log10(1-(1-e).^4);
end
    %display(sprintf('Matched features: %d \t No. of inliers: %d.\n',matched_feat, max_inliers));
   
end    

