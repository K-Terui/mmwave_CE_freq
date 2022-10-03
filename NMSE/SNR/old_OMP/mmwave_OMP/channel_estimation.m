function H_freq_est = channel_estimation(x_hat, index, K, N_r, N_t, AR_dic, AT_dic)

    % calc estimated channel (frequency domain), H_freq_est
    % size(H_freq_est) = N_r, N_t, K

    % make AoA phi_l AoD theta_l
%     Index = T_set;
    H_freq_est = zeros(N_r, N_t, K);
    for k = 0:K-1
        index_kth = index(:,k+1);
        phi = rem((index_kth-1),64)+1;
        theta = fix((index_kth-1)/64)+1;
        a_R_est = AR_dic(:,phi);
        a_T_est = AT_dic(:,theta);
    
        % calc vectorize channel
        xi_est = x_hat(index_kth);
        vec_h_est = krp(conj(a_T_est), a_R_est)*xi_est;
    
        H_freq_est(:,:,k+1) = reshape(vec_h_est, N_r, N_t);
    end


%     % calc channel gain
%     delta_vec_est = reshape(x_hat, numel(index), K);
%     
%     % calc array (kron)
%     kron_array_est = kron(eye(K), krp(conj(a_T_est),a_R_est));
%     vec_h_est = kron_array_est * reshape(delta_vec_est,K*numel(T_set),1);
%     H_freq_est = reshape(vec_h_est, N_r, N_t, K);
    
end