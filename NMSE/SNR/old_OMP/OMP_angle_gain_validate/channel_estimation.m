function H_freq_est = channel_estimation(x_hat, index, K, N_r, N_t, AR_dic, AT_dic, phi_l, theta_l, alpha_l, G_r, G_t, N_C)

    % calc estimated channel (frequency domain), H_freq_est
    % size(H_freq_est) = N_r, N_t, K

    % make AoA phi_l AoD theta_l
%     Index = T_set;
    H_freq_est = zeros(N_r, N_t, K);
    for k = 0:K-1
        index_kth = index(:,k+1);
        phi_index   = double(rem((index_kth-1),64)+1);
        theta_index = double(fix((index_kth-1)/64)+1);
        phi   = (phi_index-1)/G_r*2*pi;
        theta = (theta_index-1)/G_t*2*pi;
        a_R_est = AR_dic(:,phi_index);
        a_T_est = AT_dic(:,theta_index);
%         % 角度既知
%         phi = phi_l;
%         theta = theta_l;
%         a_R_est = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l));
%         a_T_est = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l));
%         
        % calc vectorize channel
        % 利得既知
        dft_mat = dftmtx(K); % fft(eye(K))と一緒
        p_dft = dft_mat(:,1:N_C);
        temp = (p_dft*alpha_l).';
        xi_est = temp(:,k+1);
        %         xi_est = x_hat(index_kth);
        
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