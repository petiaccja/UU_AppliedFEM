function [ c ] = fem_solver( x, a, f )
% FEM_SOLVER Simple FEM solver.
    A = a*stiffness_matrix_ddu(x, 1e7);
    b = load_vector(x, f);
    c = A\b;    
end
