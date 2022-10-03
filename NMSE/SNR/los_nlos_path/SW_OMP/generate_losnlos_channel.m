function [H_freq, PL] = generate_losnlos_channel(K, N_C, N_r, N_t, sn_dB, power_t, L, pdr)
    % write MIMO OFDM mmwave channel (frequency domain), H_freq
    % size(H_freq) = N_r, N_t, K

    PL = db2pow(power_t-sn_dB);
    
    %% decide angular parameter and array response
    % make AoA phi_l AoD theta_l
    % LOS : l = 0, NLOS : l = 1,2, ... , L-1
    % L : num. of paths (LOS + NLOS)
%     L_los   = 1;
%     L_nlos  = L - 1;
%     phi_los    = rand(1, L_los) *2*pi;
%     theta_los  = rand(1, L_los) *2*pi;
%     phi_nlos   = rand(1, L_nlos)*2*pi;
%     theta_nlos = rand(1, L_nlos)*2*pi;
    phi_l   = rand(1, L)*2*pi;
    theta_l = rand(1, L)*2*pi;

%     a_R_los  = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_los));
%     a_T_los  = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(theta_los));
%     a_R_nlos = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_nlos));
%     a_T_nlos = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_nlos));
    a_R = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l));
    a_T = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l));
    
   %% caluculate channel gain
%     alpha_los  = zeros(L_los , K);
%     alpha_nlos = zeros(L_nlos, K);
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
    
   %% generate channel
    % caluculate khatri rao product of array response matrices
    krp_array = kron(eye(K),krp(conj(a_T),a_R));
    
    % representation of channel response realizations
    vec_h = krp_array * reshape(alpha, K*L, 1);
    H_freq = reshape(vec_h, N_r, N_t, K);

end