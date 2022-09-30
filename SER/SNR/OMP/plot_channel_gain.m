% plot_channel_gain
% plot cdf of channel gain

% initialize
% clear
times = 10;             % Monte Carlo times
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
power_t = 30;             % transmitter power [dBm]

% caluculate channel gain
channel_gain = caluculate_channel_gain(times, L, N_C, N_t, N_r, K, sn_dB, power_t);

hold off
hold on

for sn = 1 : length(sn_dB) 
    
    cdf_fig = cdfplot(channel_gain(sn,:));
    cdf_fig.LineWidth = 2;
    
end

xlabel("Channel Gain (dB)")
ylabel("Comulative Distribution")
legend(["SNR = -15dB", "SNR = -10dB","SNR = -5dB", "SNR = 0dB", "SNR = 5dB", "SNR = 10dB",], "Location","southeast")
title(['CDF (K = ', num2str(K), ')'])