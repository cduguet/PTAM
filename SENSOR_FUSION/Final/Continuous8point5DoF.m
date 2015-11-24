function [omega,vel0,correl]= Continuous8point5DoF(x,y,z,u,v,w)

    %Step 1: Define matrix A constructed with the Optical Flow data
  
A = [w.*y-v.*z     u.*z-w.*x    v.*x-u.*y   x.^2    2.*x.*y    2.*x.*z   y.^2    2.*y.*z    z.^2];
 
    %Determine the SVD of A, and solve Ae=0
    [U,~,V]=svd(A);
    
    % positive  depth constraint
    if [mean(u) mean(v)]*V(1:2,end) > 0
        V=-V; U=-U;
    end
    e = V(:,end);
    e = e/norm(e(1:3));

    %recover the velocity vector and symmetric matrix  s 
    vel0 = e(1:3);
    
    %invert Z axis
    vel0(3) =-vel0(3);
    
    s = [[e(4) e(5) e(6)];
         [e(5) e(7) e(8)];
         [e(6) e(8) e(9)]];
     %     Eigenvalue decomposition of the symmetric matrix
                
    [V1,lambda] = eig(s);
    lambda = diag(lambda);

    %DEBUG: Dealing with complex eigenvalues 
    % The nearest matrix to a one with complex eingenvalues is the one with
    % real truncated igenvalues [\cite
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/290605]
    if ~isreal(lambda)
     warning('time %d: complex eigenvalues!. Truncating... \n', i);
     return
     lambda = real(lambda);
    end
    
    % eig() function does not give the eigenvalues sorted [\cite 
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/156868]. 
    % We have to do it manually
    if lambda(1) < lambda(2) || lambda(2) < lambda(3)
        [lambda,idx] = sort(lambda,1,'descend');
        V1 = V1(:,idx);
    end
    
    %in case the eigenvalues do not complain with our restrictions:
    lambda(1) = max([lambda(1) 0]);  
    lambda(3) = min([lambda(3) 0]);
    lambda(2) = max([lambda(2) lambda(3)]);
    lambda(2) = min([lambda(2) lambda(1)]);
    
    
    %Projet symmetric matrix onto the symmetric epipolar space S 
    sigma = zeros(3,1);
    sigma(1) = (2*lambda(1) + lambda(2) - lambda(3))/3;
    sigma(2) = (lambda(1) + 2*lambda(2) + lambda(3))/3;
    sigma(3) = (2*lambda(3) + lambda(2) - lambda(1))/3;
    
    %Recover velocity from the symmetric epipolar part
    Lambda = sigma(1) - sigma(3);
    Theta = acos(-sigma(2)/Lambda);
    if Theta > pi || Theta < 0
        warning('Lambda or Theta in symmetric epipolar component are out of bounds');
    end
    
    Ry  = inline(' [[cos(x) 0 sin(x)]; [0  1   0 ]; [-sin(x) 0  cos(x)]] ');
    Rz  = inline(' [[cos(x) -sin(x) 0]; [sin(x) cos(x) 0] ; [0 0 1] ]');
    V = V1 * Ry(Theta/2 - pi/2)';
    U = -V * Ry(Theta);

    SIGMA_Lambda = diag([Lambda Lambda 0]);
    SIGMA_1      = diag([1 1 0]);
    
    %construct matrix with the 4 solutions as columns
    omega = zeros(3,4);
    v_     = zeros(3,4);
    
    % Calculate the 4 different solutions:
    aux = U * Rz(pi/2) * SIGMA_Lambda * U';     omega(:,1) = [aux(8) -aux(7) aux(4)];
    aux = V * Rz(pi/2) * SIGMA_1      * V';     v_(:,1)     = [aux(8) -aux(7) aux(4)];
    
    aux = U * Rz(-pi/2) * SIGMA_Lambda * U';    omega(:,2) = [aux(8) -aux(7) aux(4)];
    aux = V * Rz(-pi/2) * SIGMA_1      * V';    v_(:,2)     = [aux(8) -aux(7) aux(4)];
    
    
    aux = V * Rz(pi/2) * SIGMA_Lambda * V';     omega(:,3) = [aux(8) -aux(7) aux(4)];
    aux = U * Rz(pi/2) * SIGMA_1      * U';     v_(:,3)     = [aux(8) -aux(7) aux(4)];
    
    
    aux = V * Rz(-pi/2) * SIGMA_Lambda * V';    omega(:,4) = [aux(8) -aux(7) aux(4)];
    aux = U * Rz(-pi/2) * SIGMA_1      * U';    v_(:,4)     = [aux(8) -aux(7) aux(4)];
    
    
    % choose the solution with v most similar to vel0
    [correl,idx] = max(vel0'*v_);
    
    omega = omega(:,idx);
%     v_ = v_(:,idx);
end