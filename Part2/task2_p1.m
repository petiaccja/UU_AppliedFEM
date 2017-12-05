% mesh
geometry = @circleg;

% exact solution
u_exact = @(x) sin(2*pi*x(1,:)).*sin(2*pi*x(2,:));

% perturbation function
f = @(x) 8*pi^2*sin(2*pi*x(1,:)).*sin(2*pi*x(2,:));


hmax = 1./(2.^(1:5))';
err_norm = zeros(length(hmax),1);

for i=1:length(hmax)
    % mesh params and generation
    [p,e,t] = initmesh(geometry, 'hmax', hmax(i));

    % assemble and solve linear system
    A = stiffness_matrix(p, e, t);
    b = load_vector(p, e, t, f);

    I = eye(length(p));
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = u_exact(p(:,e(1,:)));

    c = A\b;
    c_exact = u_exact(p)';

    % calculate error
    err = c_exact - c;
    err_norm(i) = err'*A*err;
end

% plot error
figure(1);

loglog(hmax, err_norm);

polycoeff = polyfit(log(hmax), log(err_norm), 1);
title('h vs error');
xlabel('h_{max}');
ylabel('error');
disp(['error ~ hmax^', num2str(polycoeff(1))]);

% plot solution and reference
figure(3);
trimesh(t(1:3,:)', p(1,:), p(2,:), c);
figure(4);
trimesh(t(1:3,:)', p(1,:), p(2,:), c_exact);