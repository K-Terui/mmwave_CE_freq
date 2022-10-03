function [x_hat, index] = OMP(y, A, L)

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

end