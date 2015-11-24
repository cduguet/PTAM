%Cristian Duguet

%Construct cross product matrix out of vector 


function CrossMat = VectorToCrossMat(v)
if length(v) ~=3;
    error('Vector length is not three');
end

CrossMat =      [[0     -v(3) v(2)];
                [v(3)    0   -v(1)];
                [-v(2)  v(1)    0]];
end