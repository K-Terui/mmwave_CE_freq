function [H_freq_est,gain_est] = channel_estimation(xi_tilde, K, N_r, N_t, G_t, G_r, T_set)
    % caluculate estimated channel (frequency domain), H_freq_est
    % size(H_freq_est) = N_r, N_t, K
    
    %% estimate array response matrices
    % estimate index
    Index = T_set;
    phi_index = double(rem((Index-1),G_r)+1);
    theta_index = double(fix((Index-1)/G_t)+1);
    
    % estimate AoA & AoD
    phi_est = (phi_index-1)*2*pi/G_r;
    theta_est = (theta_index-1)*2*pi/G_t;
    
    % decide array response matrices
    a_R_est = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_est'));
    a_T_est = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_est'));
    
    %% caluculate vectorized channel gain
    gain_est = reshape(xi_tilde, size(xi_tilde,1)*K,1);
    
    %% caluculate khatri rao product of array response matrices
    krp_array_est = kron(eye(K),krp(conj(a_T_est),a_R_est));
    
    %% estimate channel response realizations
    vec_h = krp_array_est * gain_est;
    H_freq_est = reshape(vec_h, N_r, N_t, K);
    
end