
function inliers_idx = EstimateInliers(x,y,u,v,H,treshold)

Hinv = inv(H);

feats_1 = [x  y  ones(size(x))]';
feats_2 = [x+u  y+v  ones(size(x))]';

diff_1 = feats_2(1:2,:) - Normalize(H*feats_1);
diff_2 = feats_1(1:2,:) - Normalize(Hinv*feats_2);

diff = sqrt(sum(diff_1.^2 + diff_2.^2,1));

inliers_idx  = find(diff <= treshold);

end

% normalize 1st and 2nd component w.r.t. the 3rd
% components rowwise
function normalized  = Normalize(points)
normalized  = points(1:2,:) ./ repmat(points(3,:),2,1);
end 
