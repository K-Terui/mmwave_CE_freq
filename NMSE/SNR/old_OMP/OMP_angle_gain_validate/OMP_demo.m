clear
K = 2;
M = 100;
N = 1000;

monte_calro = 100;

NMSE = zeros(monte_calro,1);
for itr_t = 1:monte_calro
  A = randn(M, N);
  x   = zeros(N, 1);
  x_d = randn(K, 1);
  x(randperm(N, K)) = x_d;
  
  y = A*x;
  
  x_hat = OrthogonalMatchingPursuit(y, A, K);
  
  NMSE(itr_t) = norm(x - x_hat)^2/norm(x)/N;
  plot(1:monte_calro, pow2db(NMSE))
end
