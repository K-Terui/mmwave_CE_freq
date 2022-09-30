function [xi_tilde, T_set] = MMSE(MSE_value, Upsilon_w, y_k, D_w, K, M, L_r, eta)

    %% Initialize the residual vectors to the input signal vectors
    % and support estimate
    y_w = (inv(D_w))'*y_k;
    res = y_w; % residual matrix r[k] -> res = M*L_r,1,K
    T_set = int16.empty;
    
    %% While phase start
    while MSE_value > eta
        % Distributed correlation
        c =  Upsilon_w'*res;
        % Find the maximum projection along the different spaces calc sum of correlation vector at every subcarrier
        abs_c = abs(sum(c, 3));
        [Max, index] = max(abs_c, [], 1);

        % Update the current guess of the common support
        T_set = union(T_set,index);
        
        % Project the input signal onto the subspace given by the support using WLS
%         xi_tilde = pinv(Upsilon_w(:,T_set))*y_w;
        xi_tilde = Upsilon_w(:,T_set)'*inv(Upsilon_w(:,T_set)*Upsilon_w(:,T_set)'+eye(M*L_r))*y_w;
        
        % Update residual
        res = y_w - (Upsilon_w(:,T_set)*xi_tilde);
        
        % Compute the current MSE
        MSE_value = sum(diag(res'*res))/(K*M*L_r);
        
    % End while
    end
 
end