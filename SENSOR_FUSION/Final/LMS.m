% Calculate the inliers and outliers of a group of optical Flow points
% using an Homography model

function inliers_idx = LMS(u,v,maxIter)

matched_feat = length(u);
inliers_idx = [];
if matched_feat < 9; 
    warning('not enough optical flow components, not filtering \n');

    return
end

sample_num = ceil(matched_feat/2);

min_median = Inf;

for j = 1:maxIter
            randompos = ceil(single(matched_feat)* rand(1,sample_num));
            random_u = u(randompos);
            random_v = v(randompos);
            
      %% Estimate the mean and standard deviations of the set
      maximum_ = max(random_u.^2 + random_v.^2);
      minimum_ = min(random_u.^2 + random_v.^2);
      median = maximum_ - minimum_;
%       sigma2_u = var(u);
%       sigma2_v = var(v);
        
      %% Discern between out/ inliers 
      
%       inliers_idx = find( (u-mean_u).^2 + (v-mean_v).^2  < (sigma2_u + sigma2_v));
%       inliers_no = length(inliers_idx);
      
      if median < min_median
          inliers_idx = randompos; 
          min_median = median;
          sel_max = maximum_;
          sel_min = minimum_;
      end
%       if max_inliers < inliers_no
%           max_inliers = inliers_no;
%           max_mean_u = mean_u;
%           max_mean_v = mean_v;
%           max_inliers_idx = inliers_idx;
%       end
      
%        proboutlier  = 1- max_inliers/double(matched_feat);
%        e = min(e,proboutlier);
%        loopno = log10(1-p)/log10(1-(1-e).^4);
end
quad = u.^2 + v.^2;
inliers_idx = find( (quad  < sel_max) .* (quad > sel_min));


    %display(sprintf('Matched features: %d \t No. of inliers: %d.\n',matched_feat, max_inliers));
   
end    


