function [ b ] = load_vector(p, e, t, func)
    N = size(p, 2);
    
    b = zeros(N,1);

    for i=1:size(t, 2)
       i1 = t(1, i); 
       i2 = t(2, i);
       i3 = t(3, i);
       p1 = p(:, i1);
       p2 = p(:, i2);
       p3 = p(:, i3);
       
       zi1 = func(p1) / 3;
       zi2 = func(p2) / 3;
       zi3 = func(p3) / 3;
       
       p1 = p1 - p3;
       p2 = p2 - p3;
       area = abs(p1(1)*p2(2) - p1(2)*p2(1))/2;
       
       b(i1) = b(i1) + zi1*area;
       b(i2) = b(i2) + zi2*area;
       b(i3) = b(i3) + zi3*area;
    end
end