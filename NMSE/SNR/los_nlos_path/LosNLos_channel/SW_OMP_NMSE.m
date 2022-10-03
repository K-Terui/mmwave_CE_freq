%{
SW-OMP NMSE
%}
% make received signal
clear
% definition
times = 1;             % Monte Carlo times (8/31 1280)(9/1 6400)(9/2 2000)(9/6 5600)
M = 80;                   % training flame length 
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
phi_l   = zeros(sn_length,L);
theta_l = zeros(sn_length,L);
for t = 1:times
    for sn = 1 : sn_length
        %% generate transmit pilot symbols assumed QPSK system
        % symbol_tr(Lt,1,M,K)
%         symbol_tr = generate_symbol(L_t, M, K);
        symbol_tr = pskmod(unidrnd(4,K,1)-1, 4, pi/4, 'gray');
        q_tr = pskmod(unidrnd(4,N_s,M)-1, 4, pi/4, 'gray');

        %% generate channel
        [H_freq, phi_l_sn, theta_l_sn, gain_sn, PL] = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);
        % genie aided generate channel
%         [H_freq, phi_l, theta_l, alpha, PL, T_set_genie] = genie_generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L, G_r, G_t);

        

    end
end

