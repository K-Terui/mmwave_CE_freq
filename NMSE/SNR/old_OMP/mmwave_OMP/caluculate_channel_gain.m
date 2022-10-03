function channel_gain = caluculate_channel_gain(times, L, N_C, N_t, N_r, K, sn_dB, power_t)
% to caluculate channel gain of frequency selective mm-wave channel (freqency domain)
% Input : times is Monte Carlo times, L is num of paths, N_C is delay taps
%         N_t is num of transmitter antennas, N_r is num of receiver antennas
%         K is OFDM subcarriers 16, sn_dB is SNR roop iterate (-15,-10,-5,0,5,10)[dB]
%         power_t is transmitter power [dBm]
%
% Output : channel_gain is matrix that contains channel gain. (length(sn),K*times)

    % initialize
    store_fro_true = zeros(times,K);
    channel_gain = zeros(length(sn_dB),times*K);
    
    for sn = 1 : length(sn_dB)
        
        for t = 1:times
        % generate channle

                H_freq = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);

                % caluculate channel gain
                true_vec = reshape(H_freq, N_r*N_t, K);
                fro_true = vecnorm(true_vec).^2;
                store_fro_true(t,:) = pow2db(fro_true);
                
        end
    
        channel_gain(sn,:) = reshape(store_fro_true,1,K*times);
        
    end

end