% monte-carlo integrates phii*phij

N = 10000000;
points = rand(2, N);

f1 = @(x) [-1, -1]*x + 1;
f2 = @(x) [1, 0]*x + 0;
f3 = @(x) [0, 1]*x + 0;

for i=1:N
   p = points(:, i);
   % reflect points to keep it in unit triangle instead of unit square
   if (p(1) + p(2) > 1)
       points(:, i) = [1; 1] - p;
   end
end

M = zeros(3,3);

for i=1:3 
   for j=i:3
      if (i==1)
          a = f1;
      elseif (i==2)
          a = f2;
      else
          a = f3;
      end
      if (j==1)
          b = f1;
      elseif (j==2)
          b = f2;
      else
          b = f3;
      end
      f = @(x) a(x).*b(x);
      r = f(points);
      int = mean(r);
      M(i,j) = int;
      M(j,i) = int;
   end
end
%scatter(points(1,:),points(2,:));
disp(M);