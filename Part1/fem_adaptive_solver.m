function [ x, c, e, r, err ] = fem_adaptive_solver( x, a, f, tolerance, maxPoints, e )
% FEM_ADAPTIVE_SOLVER Iteratively refines the mesh to reach a predefined a
%   posteriori error estimate.
%   Returns: vertices, coefficients, estimated total error vs iteration,
%   residual, error contribution of points

    N = length(x)-1;
    
    % calculate fem solution
    A = a*stiffness_matrix_ddu(x, 1e7);
    b = load_vector(x, f);
    c = A\b;
    
    % calculate discrete laplacian
    M = mass_matrix(x);
    p = (-A*c);
    ddu = (a*M)\p;
    
    % make a posteriori error estimate
    err = zeros(N, 1);
    for i = 1:N
        h = x(i+1) - x(i);
        lerpddu = @(x0) ((x0-x(i))*ddu(i+1) + (x(i+1)-x0)*ddu(i)) / (x(i+1)-x(i)); 
        err(i) = h^2*integral(@(x0) (f(x0) + a*lerpddu(x0)).^2, x(i), x(i+1));
    end
    
    % residual
    r = a*ddu + f(x)';
    
    % calculate error and collect it for each iteration
    sumerr = sum(err);
    e = [e, sumerr];
    
    % refine mesh and do again
    if (sumerr > tolerance && length(x) < maxPoints) 
        xref = refine_mesh1(x, err);
        [x, c, e, r, err] = fem_adaptive_solver(xref, a, f, tolerance, maxPoints, e);
    end    
end


function [M] = mass_matrix(x)
    N = length(x)-1;
    
    M = zeros(N+1, N+1);
    for i = 1:N
        hi = x(i+1) - x(i);
        ii = [i, i+1];
        M(ii,ii) = M(ii,ii) + [1/3, 1/6; 1/6, 1/3]*hi;
    end
end


% my mesh refining function cause the one in the description seemed too
% trivial
function [xref] = refine_mesh1(x, err) 
    N = length(x)-1;
    m = mean(err);
    splitcount = floor(sqrt(err/m) + 1);
    
    xref = [];
    for i=1:N
        a = x(i);
        b = x(i+1);
        h = (b-a)/splitcount(i);
        xref = [xref, a:h:(b-h/2)];
    end
    xref = [xref, x(N+1)];
end


% the mesh refining function given in the project description
function [xref] = refine_mesh2(x, err) 
    lambda = 0.9;
    m = max(err);
    for i = 1:length( err )
        if err(i) > lambda*m
            x = [x (x( i +1)+ x( i ))/2];
        end
    end
    xref = sort(x);
end

