function [ deltaM ] = mass_loss( p, t, c_initial, c )
%MASS_LOSS Summary of this function goes here
%   Detailed explanation goes here
    diff = c_initial - c;
    
    sum = 0;
    
    for i=1:size(t, 2)
       d = diff(t(1:3, i));
       
       i1 = t(1, i); 
       i2 = t(2, i);
       i3 = t(3, i);
       p1 = p(:, i1);
       p2 = p(:, i2);
       p3 = p(:, i3);       
       p1 = p1 - p3;
       p2 = p2 - p3;
       area = abs(p1(1)*p2(2) - p1(2)*p2(1))/2;
       
       sum = sum + mean(d)*area;        
    end

    deltaM = sum;
end

