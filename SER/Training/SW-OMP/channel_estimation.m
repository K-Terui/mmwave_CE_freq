function H_freq_est = channel_estimation(xi_tilde, T_set, K, N_r, N_t, AR_dic, AT_dic, G_r, G_t, phi_l, theta_l)

    % calc estimated channel (frequency domain), H_freq_est
    % size(H_freq_est) = N_r, N_t, K

    % make AoA phi_l AoD theta_l
    index = T_set;
    phi_index   = double(rem((index-1),64)+1);
    theta_index = double(fix((index-1)/64)+1);
    phi   = (phi_index-1)/G_r*2*pi;
    theta = (theta_index-1)/G_t*2*pi;
%     % 角度既知
%     phi = phi_l(:,1:numel(T_set));
%     theta = theta_l(:,1:numel(T_set));
    
%     a_R_est = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi));
%     a_T_est = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta));
    a_R_est = sqrt(1/N_r)*exp(1i*pi*(0:N_r-1)'.*cos(phi));
    a_T_est = sqrt(1/N_t)*exp(1i*pi*(0:N_t-1)'.*cos(theta));
%     a_R_est = AR_dic(:,phi_index);ra
%     a_T_est = AT_dic(:,theta_index);

    % calc channel gain
    delta_vec_est = reshape(xi_tilde, numel(T_set), K);
    
    % calc array (kron)
    kron_array_est = kron(eye(K), krp(conj(a_T_est),a_R_est));
    vec_h_est = kron_array_est * reshape(delta_vec_est,K*numel(T_set),1);
    H_freq_est = reshape(vec_h_est, N_r, N_t, K);
    
end