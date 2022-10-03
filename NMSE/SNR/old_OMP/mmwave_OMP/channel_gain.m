%% 
% initialize
clear
times = 100;             % Monte Carlo times
M = 80;                    % training flame length 
N_s = 1;                  % data stream Ns=Lt(p2949 leftside)
sample = 1000;            % num of sample
sps = 4;                  % symbol per sample
L = 4;                    % num of paths 
N_C = 4;                  % delay taps
Ts = 1/1760*10^(-6);      % sampling period
K = 16;                   % OFDM subcarriers 16
N_t = 32;                 % num of transmitter antennas
N_r = 32;                 % num of receiver antennas
L_t = 1;                  % num of transmitter RF chain
L_r = 4;                  % num of receiver RF chain
G_t = 64;                 % grid of size for the AoA
G_r = 64;                 % grid of size for the AoD
freq = 60*10^(9);         % mm-wave freq (802.11ad says 60Ghz)
beta_r = 0.8;             % rolloff factor of RCF
var = 1;                  % variance of noise
sn_dB = -15:5:10;         % SNR roop iterate (-15,-10,-5,0,5,10)[dB]
power_t = 30;             % transmitter power [dBm]
s_pre = rng;              % pseudorandomly precoder
s_com = rng;              % pseudirabdinly combiner
N_Q = 2;                  % quantization bits
store_fro_true = zeros(times,K);
cdf_channle_gain = zeros(length(sn_dB),times*K);

for sn = 1 : length(sn_dB)
    for t = 1:times
    %%
    % generate channle

            H_freq = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);

            %%
            % channel gain
            true_vec = reshape(H_freq, N_r*N_t, K);
            fro_true = vecnorm(true_vec).^2;
            store_fro_true(t,:) = pow2db(fro_true);
    end
    cdf_channle_gain(sn,:) = reshape(store_fro_true,1,K*times);

end

hold off
hold on

for sn = 1 : length(sn_dB) 
    cdf_fig = cdfplot(cdf_channle_gain(sn,:));
    cdf_fig.LineWidth = 2;
end
xlabel("Channel Gain (dB)")
ylabel("Comulative Distribution")
legend(["SNR = -15dB", "SNR = -10dB","SNR = -5dB", "SNR = 0dB", "SNR = 5dB", "SNR = 10dB",], "Location","southeast")
title(['CDF (M = ', num2str(M), ')'])