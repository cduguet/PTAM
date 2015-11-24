function normalized=quat_normalize(quaternion)
if length(quaternion) ~= 4  
    error('Quaternion does not have the format [rx ry rz w]!!');
    return;
end
normalized = quaternion/norm(quaternion);
end
