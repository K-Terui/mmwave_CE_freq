function [xi_tilde, T_set] = genie_SW_OMP(y_k, Upsilon_w, D_w, T_set_genie)
    
    %% Initialize the residual vectors to the input signal vectors
    % and support estimate
    y_w = (inv(D_w))'*y_k;
%     res = y_w; % residual matrix r[k] -> res = M*L_r,1,K
    T_set = T_set_genie;

    %% only gain estimation
    % Project the input signal onto the subspace 
    % given by the support using WLS
    xi_tilde = pinv(Upsilon_w(:,T_set))*y_w;

%     % Update residual
%     res = y_w - (Upsilon_w(:,T_set)*xi_tilde);
%     
%     % Compute the current MSE
%     MSE_value = sum(diag(res'*res))/(K*M*L_r);

end