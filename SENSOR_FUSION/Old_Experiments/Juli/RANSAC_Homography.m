% Calculate the inliers and outliers of a group of optical Flow points
% using an Homography model

function max_inliers_idx = RANSAC_Homography(x,y,u,v,maxIter,sample_num,RANSAC_treshold)

max_inliers =0;
matched_feat = length(x);

for j = 1:maxIter
    
    
        coflag=1;
        while coflag          
            randompos = ceil(single(matched_feat)* rand(1,sample_num));
            random_x1 = x(randompos);
            random_y1 = y(randompos);
            random_x2 = x(randompos) + u(randompos);
            random_y2 = y(randompos) + v(randompos);
            
            if CheckColinearity([random_x1; random_y1; ones(size(random_x1))]) ==0
                coflag=0;
            end
        end
    
      %%MAKE HOMOGRAPHY MATRIX
      B = zeros(2*sample_num,9);
      B(1:2:2*sample_num,1:3) = [random_x1 random_y1 ones(size(random_x1))]; 
      B(2:2:2*sample_num,4:6) = [random_x1 random_y1 ones(size(random_x1))];
      B(1:2:2*sample_num,7:9) = -[random_x1 random_y1 ones(size(random_x1))] .* repmat(random_x2,1,3);
      B(2:2:2*sample_num,7:9) = -[random_x1 random_y1 ones(size(random_x1))] .* repmat(random_y2,1,3);
      
      [~,~,VV] = svd(B);
      
      H = [VV(1:3,9) VV(4:6,9) VV(7:9,9)]';
      
      % -----------------------------------------------------Cluttered code
      %       fprintf('H before normalization: cond = %d, norm = %d\t', cond(H), norm(H));
      %       H = T2\H*T1;
      %       fprintf('H after normalization: cond = %d, norm = %d\n', cond(H), norm(H));
      % H shows a worst condition after normalization.
      %--------------------------------------------------------------------

      inliers_idx = EstimateInliers(x,y,u,v,H,RANSAC_treshold);
      inliers_no = length(inliers_idx);
      if max_inliers < inliers_no
          max_inliers = inliers_no;
          max_H = H;
          max_inliers_idx = inliers_idx;
      end
%        proboutlier  = 1- max_inliers/double(matched_feat);
%        e = min(e,proboutlier);
%        loopno = log10(1-p)/log10(1-(1-e).^4);
    end
    
    %display(sprintf('Matched features: %d \t No. of inliers: %d.\n',matched_feat, max_inliers));

end



% % Estimate which points of  (x,y) ->(x+u,y+v) can be inliers of the
% transformation defined by the Homography matrix H 

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

%Check colinearity of columns
function lol=CheckColinearity(Mat)
no_cols = size(Mat,2);
for i=1:no_cols-2
    for j=i+1:no_cols-1
        for k=j+1:no_cols
            if Mat(:,i)'*cross(Mat(:,j),Mat(:,k)) ==0
                lol=1;
                return
            end
        end
    end
end
lol=0;
return
end




    
    
    