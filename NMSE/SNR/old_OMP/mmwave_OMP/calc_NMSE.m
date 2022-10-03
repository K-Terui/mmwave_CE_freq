function NMSE_dB = calc_NMSE(H_freq_est, H_freq, K, N_r, N_t)

    % calc Frobenius norm
    error = H_freq_est - H_freq;
    error_vec = reshape(error, N_r*N_t, K);
    fro_error = vecnorm(error_vec).^2;
    true_vec = reshape(H_freq, N_r*N_t, K);
    fro_true = vecnorm(true_vec).^2;
    
    % calc NMSE
    numerator = sum(fro_error);
    denominator = sum(fro_true);
    NMSE_dB = pow2db(numerator/denominator);

end