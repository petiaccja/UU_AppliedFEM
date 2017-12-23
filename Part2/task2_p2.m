hmax_v = [1/5, 1/20];
for rep=1:2
    
% mesh
geometry = @circleg;
hmax = hmax_v(rep);
[p,e,t] = initmesh(geometry, 'hmax', hmax);
N = length(p);

% properties
f = @(x) 0;
alpha = 0.01;
R = 0.5;
r = 0.3;
deltaT = 0.05;
Tend = 30;
rho = 10;

% initial condition
c_initial = zeros(N, 1);
for i=1:N
   pt = p(:, i);
   if (abs(R - norm(pt)) < r) 
      c_initial(i) = rho;
   end
end


% assemble and solve linear system
M = mass_matrix(p, e, t);
A = alpha*stiffness_matrix(p, e, t);
b = load_vector(p, e, t, f);

G = (0.5/deltaT)*M + (1/2)*A;
H = (-0.5/deltaT)*M + (1/2)*A;

I = speye(N);
G(e(1,:),:) = I(e(1,:),:);

% plot initial state
figure(1);
trimesh(t(1:3,:)', p(1,:), p(2,:), c_initial);
title('Initial conditions');

% iterate solver
c_prev = c_initial;
figure(2);
animplot = trimesh(t(1:3,:)', p(1,:), p(2,:), c_initial);
ml = zeros(1);

iter = 1;
for time=0:deltaT:Tend
    % solve system
    d = b - H*c_prev;   
    d(e(1,:)) = 0;
    
    %c = G\d;
    [c, flag] = bicg(G, d);
    c_prev = c;
    
    % calculate mass loss
    ml(iter) = mass_loss( p, t, c_initial, c);
    
    % plot animation
	if (mod(iter, 10) == 1 && rep==2)
        trimesh(t(1:3,:)', p(1,:), p(2,:), c);
        axis manual
        axis([-1, 1, -1, 1, min(c_initial)-1, max(c_initial)+1]);  
        title(['T = ', num2str(time)]);
        drawnow;
    end
    
    iter = iter + 1;
end
set(gcf,'visible','off');


if (rep == 2)
figure(3);
trimesh(t(1:3,:)', p(1,:), p(2,:), c_prev);
title('Final solution');
end 

figure(3+rep);
plot(0:deltaT:Tend, ml);
restext = [' (low res) '; ' (high res)'];
title(['Mass loss vs time', restext(rep, :)]);
xlabel('Time');
ylabel('Mass loss');


end


