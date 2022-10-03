function [H_freq, phi_l, theta_l, alpha, PL, T_set_genie] = genie_generate_losnlos_channel(K, N_C, N_r, N_t, sn_dB, power_t, L, pdr, G_r, G_t)
    % write MIMO OFDM mmwave channel (frequency domain), H_freq
    % size(H_freq) = N_r, N_t, K
    % genei aided (known parameter: phi_l, theta_l, alpha)

    PL = db2pow(power_t-sn_dB);
    
    %% decide angular parameter and array response
    % make AoA phi_l AoD theta_l
    % genie aided
    phi_index   = randi(G_r, 1, L);
    theta_index = randi(G_t, 1, L);
    phi_l   = (phi_index-1)*2*pi/G_r;
    theta_l = (theta_index-1)*2*pi/G_r;
    T_set_genie = theta_index*G_t-G_t+1+(phi_index-1);

    a_R = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l));
    a_T = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l));
    
    %% caluculate channel gain
    alpha_losnlos = zeros(L, K);
    dft_mat = dftmtx(K)/sqrt(K);
    p_dft = dft_mat(:, 1:N_C);
    pow_dist = 10^(-pdr);
    nlos_pow = repmat(pow_dist, L-1, 1);
    pow_distm = cat(1,eye(1),nlos_pow);
    
    for l = 1:L
        cov_mat = p_dft*eye(N_C)/PL/L*p_dft';
        alpha_losnlos(l,:) = sqrtm(cov_mat)*crandn(K, 1);
    end
    alpha = alpha_losnlos.*pow_distm;   
    
    % caluculate khatri rao product of array response matrices
    krp_array = kron(eye(K),krp(conj(a_T),a_R));
    
    % representation of channel response realizations
    vec_h = krp_array * reshape(alpha, K*L, 1);
    H_freq = reshape(vec_h, N_r, N_t, K);

end