% plot_channel_gain
% plot cdf of channel gain

%% initialize
% clear
times = 1280;             % Monte Carlo times
L = 4;                    % num of paths 
N_C = 4;                  % delay taps
K = 16;                   % OFDM subcarriers 16
N_t = 32;                 % num of transmitter antennas
N_r = 32;                 % num of receiver antennas
L_t = 1;                  % num of transmitter RF chain
L_r = 4;                  % num of receiver RF chain
G_t = 64;                 % grid of size for the AoA
G_r = 64;                 % grid of size for the AoD
sn_dB = -15:5:10;         % SNR roop iterate (-15,-10,-5,0,5,10)[dB]
sn_length = length(sn_dB);% SNR roop length 
power_t = 30;             % transmitter power [dBm]

store_fro_true = zeros(sn_length, K, times);

%% caluculate channel gain
parfor t = 1 : times
    for sn = 1 : sn_length
        
        H_freq = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);

        % caluculate channel gain
        true_vec = reshape(H_freq, N_r*N_t, K);
        fro_true = vecnorm(true_vec).^2;
        store_fro_true(sn,:,t) = pow2db(fro_true);
        
    end
end

channel_gain= reshape(store_fro_true, sn_length, K*times);


%% plot 
hold off
figure
hold on

for sn = 1 : length(sn_dB) 
    
    cdf_fig = cdfplot(channel_gain(sn,:));
    cdf_fig.LineWidth = 2;
    
end

xlabel("Channel Gain (dB)")
ylabel("Comulative Distribution")
legend(["SNR = -15dB", "SNR = -10dB","SNR = -5dB", "SNR = 0dB", "SNR = 5dB", "SNR = 10dB",], "Location","southeast")
% title(['CDF (K = ', num2str(K), ')'])