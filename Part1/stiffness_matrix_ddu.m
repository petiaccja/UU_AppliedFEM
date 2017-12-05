function [A] = stiffness_matrix_ddu(x, boundaryHuge)
% stiffness_matrix_ddu calculates the stiffnedd matrix for -u''=f
    N = length(x)-1;
    
    A = zeros(N+1, N+1);
    for i = 1:N
        hi = x(i+1) - x(i);
        ii = [i, i+1];
        A(ii,ii) = A(ii,ii) + [1, -1; -1, 1]/hi;
    end
    A(1,1) = boundaryHuge;
    A(N+1, N+1) = boundaryHuge;
end