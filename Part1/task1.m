% Clear figures
for i=1:5
   figure(i);
   clf(i);
end

% Problem definition
x = -1:0.2:1;
f1 = @(x) 1+x;
a = 0.01;

% Solve with original nodes
figure(1);
c = fem_solver(x, a, @f2);

% Solve adaptively
[xf, cf, ef, rf, errdistf] = fem_adaptive_solver(x, a, @f2, 1e-4, 1e4, []);

% Plot previous two solutions
figure(1);
%plot(x, c);
hold on;
plot(xf, cf);
title('Solution');
xlabel('x');
ylabel('u(x)');
%legend('Simple', 'Adaptive');

% Plot the error estimate of the adaptive solver versus iteration count
figure(4);
plot(log10(ef));
title('Error estimate');
xlabel('Adaptive solver iteration count');
ylabel('Error estimate (log scale)');

figure(2);
plot((xf(1:(length(xf)-1))+xf(2:length(xf)))/2, sqrt(errdistf));
title('Error contribution');
xlabel('Element position');
ylabel('Element error contribution');

figure(3);
plot((xf(1:(length(xf)-1))+xf(2:length(xf)))/2, 1./diff(xf));
title('Mesh density');
xlabel('Element position');
ylabel('Mesh density');

figure(5);
plot(xf, rf);
xlabel('x');
ylabel('f(x) + a*u''''(x)');
title('Residual');

% another perturbing function for pde
function y = f2(x)
    R = 0.5;
    p = 10;
    r = 0.3;
    y = x;
    for i=1:length(x)
        if (abs(R - abs(x(i))) < r)
            y(i) = p;
        else
            y(i) = 0;
        end
    end
end
