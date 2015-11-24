% Calculate the inliers and outliers of a group of optical Flow points
% using an Homography model

function inliers_idx = HistographicRejection(u,v,maxIter)

matched_feat = length(u);
inliers_idx = [];
if matched_feat < 9; 
    warning('not enough optical flow components, not filtering \n');

    return
end

sample_num = ceil(matched_feat/2);

% min_median = Inf;

 maxu = max(u);
 minu = min(u);
 edgesu = linspace(minu,maxu,6);
 [nu,binu] = histc(u',edgesu);
 [~,modeu] = max(nu);
 
 try inliersu = (u<edgesu(modeu+1)) .*( u>=edgesu(modeu)); 
 catch
    fprintf('weird OF...\n');
    inliers_idx = 1:length(u);
    return
 end
 maxv = max(v);
 minv = min(v);
 edgesv = linspace(minv,maxv,6);
 [nv,binv] = histc(v',edgesv);
 [~,modev] = max(nv);
 
 inliersv = (v<edgesv(modev+1)) .*( v>=edgesv(modev));
 
 inliers_idx = find(inliersu .*inliersv);
   
end    


