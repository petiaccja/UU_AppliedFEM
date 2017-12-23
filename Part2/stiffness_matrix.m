function [ A ] = stiffness_matrix( p, e, t )
%STIFFNESS_MATRIX Summary of this function goes here
%   Detailed explanation goes here
    N = size(p, 2);
    
    A = sparse(N, N);
    
    % integrate gradient for elements
    for i=1:size(t, 2)
       i1 = t(1, i); 
       i2 = t(2, i);
       i3 = t(3, i);
       p1 = p(:, i1);
       p2 = p(:, i2);
       p3 = p(:, i3);
       v11 = integrate_dphii_dphij(p1, p2, p3, 1, 1);
       v12 = integrate_dphii_dphij(p1, p2, p3, 1, 2);
       v13 = integrate_dphii_dphij(p1, p2, p3, 1, 3);
       v22 = integrate_dphii_dphij(p1, p2, p3, 2, 2);
       v23 = integrate_dphii_dphij(p1, p2, p3, 2, 3);
       v33 = integrate_dphii_dphij(p1, p2, p3, 3, 3);
              
       A(i1, i2) = A(i1, i2) + v12;
       A(i1, i3) = A(i1, i3) + v13;
       A(i2, i3) = A(i2, i3) + v23;
       
       A(i1, i1) = A(i1, i1) + v11;
       A(i2, i2) = A(i2, i2) + v22;
       A(i3, i3) = A(i3, i3) + v33;
       
       A(i2, i1) = A(i2, i1) + v12;
       A(i3, i1) = A(i3, i1) + v13;
       A(i3, i2) = A(i3, i2) + v23;
    end
end


function [v] = integrate_dphii_dphij(p1, p2, p3, i, j)
    M = get_element_transform(p1, p2, p3);
    MiT = inv(M)';
    scale = det(M);
    
    dphi = [-1, 1, 0;...
            -1, 0, 1];
        
    d = (MiT*dphi(:,i))' * (MiT*dphi(:,j)) / 2;
    v = scale*d;
end

