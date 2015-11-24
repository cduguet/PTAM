
% Multiplication of one or many quatiernions ordered ina  row. 
% Order x y z w
function Qprod=quatmult(Q0,Q1)
% Q0 and Q1 shoud be in this order
% Q0=[x0 y0 z0 w0] % w0 is scalar part, x0,y0,z0 are vector part
% Q1=[x1 y1 z1 w1] % w1 is scalar part, x1,y1,z1 are vector part
% Multiplication is not commutative in that the products Q0Q1 and Q1Q0 are
% not necessarily equal.
N=size(Q0,1);
if  N~= size(Q1,1)
   error('Dimesions do not match!');
end

Qprod = zeros(N,4);

for i=1:N
    w0=Q0(i,4); x0=Q0(i,1); y0=Q0(i,2); z0=Q0(i,3); 
    w1=Q1(i,4); x1=Q1(i,1); y1=Q1(i,2); z1=Q1(i,3); 


%    
%     xr=(w0.*x1 + x0.*w1 - y0.*z1 + z0.*y1);
%     yr=(w0.*y1 - x0.*z1 + y0.*w1 + z0.*x1);
%     zr=(w0.*z1 + x0.*y1 - y0.*x1 + z0.*w1);
%     wr=(w0.*w1 - x0.*x1 - y0.*y1 - z0.*z1);
%     
% %     
%     According to standard from MATLAB guide 
    xr= w0.*x1 + z0.*y1 - y0.*z1 + x0.*w1;
    yr=-z0.*x1 + w0.*y1 + x0.*z1 + y0.*w1;
    zr= y0.*x1 - x0.*y1 + w0.*z1 + z0.*w1;
    wr=-x0.*x1 - y0.*y1 - z0.*z1 + w0.*w1;
    
    
%      q4.*p1 + q3.*p2 -q2.*p3 + q1.*p4;
%     -q3.*p1 + q4.*p2 + q1.*p3 + q2.*p4;
%      q2.*p1 - q1.*p2 +q4.*p3 + q3.*p4;
%     -q1.*p1 - q2.*p2 -q3.*p3 + q4.*p4;
%     
    Qprod(i,:)=[xr yr zr wr];  
end

return 
 % wr is scalar part, xr, yr, zr are vector part

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % modified by: Cristian Duguet
% % % --------------------------------

% Technical Reference: Ken Shoemake, "Animating Rotations with Quaternion Curves"
