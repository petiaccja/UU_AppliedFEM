function [b] = load_vector(x, f)
% load_vector calculates the load vector for g(u) = f(x)
    N = length(x)-1;
    
    b = zeros(N+1, 1);    
    for i = 1:N
       hi = x(i+1) - x(i);
       ii = [i, i+1];
       b(ii) = b(ii) + [f(x(i)); f(x(i+1))]*hi/2;
    end
end