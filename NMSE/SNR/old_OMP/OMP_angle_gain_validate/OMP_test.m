clear
L = 2;
M = 100;
N = 1000;

monte_calro = 100;

NMSE = zeros(monte_calro,1);
for itr_t = 1:monte_calro
    A = randn(M, N);
    x   = zeros(N, 1);
    x_d = randn(L, 1);
    x(randperm(N, L)) = x_d;

    y= A*x;
    x_hat = zeros(size(A,2),1);
    index = zeros(L, 1);
    A_hat = zeros(size(A));
    r = y;

    for i = 1:L
        % decide index
        [~, index(i)] = max(abs(A'*r));
    
        % obtain column vector from measurement matrix
        A_hat(:, i) = A(:, index(i));
    
        % clean up matrix A
        A(:, index(i)) = 0;
    
        % LS
        x_tilde = A_hat(:, 1:i)\y;
    
        % apdate residual
        r = y - A_hat(:, 1:i)*x_tilde;
    end
    % reconstruct sparse signal
    x_hat(index) = x_tilde;
    
    NMSE(itr_t) = norm(x - x_hat)^2/norm(x)/N;
    plot(1:monte_calro, pow2db(NMSE))
end