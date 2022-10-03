function [x_hat, x_tilde] = OrthogonalMatchingPursuit(y, A, L)
    %Initialize
    A_hat   = zeros(size(A));
    x_hat   = zeros(size(A,2),1);
    x_tilde = zeros(L,1);
    r = y;
    
        for i = 1:L
            [~, x_tilde(i)]   = max(abs(A.'*r));
            A_hat(:, i)       = A(:, x_tilde(i));
            A(:, x_tilde(i))  = 0;
            
            estimate_vec = A_hat(:, 1:i)\y;
            r = y - A_hat(:, 1:i)*estimate_vec;
        end
    
    x_hat(x_tilde) = estimate_vec;
end