n = 20;
A = 500;
a = zeros(n);

parfor i = 1:n
  a(i) = max(abs(eig(rand(A))));
end
