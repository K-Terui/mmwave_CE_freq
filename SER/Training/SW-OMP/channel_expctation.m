% make received signal
clear
% definition
times = 1000;             % Monte Carlo times (8/31 1280)(9/1 6400)
M = 100;                    % training flame length 
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
variance = 1;             % variance of noise
sn_dB = -15:5:10;         % SNR roop iterate (-15,-10,-5,0,5,10)[dB]
% sn_dB = 0;                % SNR roop iterate (-15,-10,-5,0,5,10)[dB]
sn_length = length(sn_dB);% SNR roop length 
power_t = 30;             % transmitter power [dBm]
s_pre = rng;              % pseudorandomly precoder
s_com = rng;              % pseudirabdinly combiner
N_Q = 2;                  % quantization bits

%initialize
norm_store = zeros(times,1);

sn_dB = sn_dB(4);


for t = 1:times
    [H_freq, phi_l_sn, theta_l_sn, gain_sn, PL] = generate_channel(K, N_C, N_r, N_t, sn_dB, power_t, L);
    true_vec = reshape(H_freq, N_r*N_t, K);
    fro_true = vecnorm(true_vec).^2;
    norm = sum(fro_true);
    norm_store(t,:) = norm;
end

expt_norm = sum(norm_store)/times;



