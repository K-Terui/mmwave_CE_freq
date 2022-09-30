function [H_freq_est_3,gain_est] = channel_estimation_3(xi_tilde, K, N_r, N_t, phi_l_sn, theta_l_sn, G_t, G_r, T_set)
    % calc estimated channel (frequency domain), H_freq_est
    % size(H_freq_est) = N_r, N_t, K
    
    %% new approach using Heath mmwave wideband eq(28) (8/28)
    % estimate array response matrices
    % estimate index
    Index = T_set;
    phi_index = double(rem((Index-1),64)+1);
    theta_index = double(fix((Index-1)/64)+1);
    
    % estimate AoA & AoD
    phi_est = (phi_index-1)*2*pi/G_r;
    theta_est = (theta_index-1)*2*pi/G_t;
    
    % decide array response matrices
    a_R_est = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_est'));
    a_T_est = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_est'));
    a_R = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l_sn));
    a_T = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l_sn));
    
    % caluculate vectorized channel gain
    gain_est = reshape(xi_tilde, size(xi_tilde,1)*K,1);
    
    % caluculate khatri rao product of array response matrices
%     krp_array_est = kron(eye(K),krp(conj(a_T_est),a_R_est));
    krp_array_est = kron(eye(K),krp(conj(a_T),a_R));    

    % representation of channel response realizations
%     vec_h = krp_array_est * gain_est;
    vec_h = krp_array_est * gain_est;
    H_freq_est_3 = reshape(vec_h, N_r, N_t, K);
    

    
end