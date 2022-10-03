function [H_freq, phi_l, theta_l, gain] = generate_channel(K, N_C, N_r, N_t, sn_dB, power_t, L)
    % write MIMO OFDM mmwave channel (frequency domain), H_freq
    % size(H_freq) = N_r, N_t, K

    PL = db2pow(power_t-sn_dB);

    % make AoA phi_l AoD theta_l
    phi_l = rand(1,L)*2*pi;
    theta_l = rand(1,L)*2*pi;
    a_R = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l));
    a_T = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l));
    
    %% new approach using Heath mmwave wideband eq(28) (8/28)
    % caluculate vectorized channel gain
    alpha = crandn(L,1)/sqrt(PL*L);
    beta = repmat(sum(exp(-1i*2*pi*(0:(K-1)).'*(0:(N_C-1))/K).'), L ,1);
    gain = alpha.*beta;
    
    % caluculate khatri rao product of array response matrices
    krp_array = krp(conj(a_T),a_R);
    
    % representation of channel response realizations
    vec_h = krp_array * gain;
    H_freq = reshape(vec_h, N_r, N_t, K);

    %% conventional approach
%     % calc channel gain
% %     alpha_l = repmat(crandn(1,L)/sqrt(PL*L), [N_C,1]);
%     alpha_l = crandn(N_C, L)/sqrt(PL*L);
% 
%     dft_mat = dftmtx(K); % fft(eye(K))と一緒
%     p_dft = dft_mat(:,1:N_C);
%     delta_vec = (p_dft*alpha_l).';    
%     
%     % calc array (kron)
%     kron_array = kron(eye(K), krp(conj(a_T),a_R));
%     vec_h = kron_array * reshape(delta_vec,K*L,1);
%     H_freq_conv = reshape(vec_h, N_r, N_t, K);

end