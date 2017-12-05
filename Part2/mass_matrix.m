function [ M ] = mass_matrix( p, e, t )
%MASS_MATRIX Summary of this function goes here
%   Detailed explanation goes here
    N = size(p, 2);
    
    M = sparse(N, N);
    
    % integrate gradient for elements
    for i=1:size(t, 2)
       i1 = t(1, i); 
       i2 = t(2, i);
       i3 = t(3, i);
       p1 = p(:, i1);
       p2 = p(:, i2);
       p3 = p(:, i3);
       v11 = integrate_phii_phij(p1, p2, p3, 1, 1);
       v12 = integrate_phii_phij(p1, p2, p3, 1, 2);
       v13 = integrate_phii_phij(p1, p2, p3, 1, 3);
       v22 = integrate_phii_phij(p1, p2, p3, 2, 2);
       v23 = integrate_phii_phij(p1, p2, p3, 2, 3);
       v33 = integrate_phii_phij(p1, p2, p3, 3, 3);
              
       M(i1, i2) = M(i1, i2) + v12;
       M(i1, i3) = M(i1, i3) + v13;
       M(i2, i3) = M(i2, i3) + v23;
       
       M(i1, i1) = M(i1, i1) + v11;
       M(i2, i2) = M(i2, i2) + v22;
       M(i3, i3) = M(i3, i3) + v33;
       
       M(i2, i1) = M(i2, i1) + v12;
       M(i3, i1) = M(i3, i1) + v13;
       M(i3, i2) = M(i3, i2) + v23;
    end
end


function [v] = integrate_phii_phij(p1, p2, p3, i, j)
    M = get_element_transform(p1, p2, p3);
    scale = det(M);
    
    Ipre = [...
        1/6    1/12    1/12;...
        1/12    1/6    1/12;...
        1/12    1/12    1/6;...
    ];
        
    d = Ipre(i,j);
    v = scale*d;
end
