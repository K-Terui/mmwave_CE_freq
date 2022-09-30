function SER = caluculate_SER(K, sn_dB, N_t, N_r, H_freq_est)

%     log_part = 1 + sn_dB*W_tr'*H_freq_est*F_tr*F_tr'*H_freq_est'*W_tr;
    % initialize
    S = zeros(N_r, N_t);
    eigen = zeros(1,K);
    
    for k = 1:K
        [~, S(:,:,k), ~] = svd(H_freq_est(:,:,k));
        eigen(1,k) = S(1,1,k);
    end
    
    SER = sum(log2(1+(db2pow(sn_dB))*eigen.^2))/K;
    
end