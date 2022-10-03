%{
SW-OMP NMSE
%}
% make received signal
clear
% definition
times = 11520;             % Monte Carlo times (8/31 1280)(9/1 6400)(9/2 2000)(9/6 5600)
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
pdr = 2.3;                % power distribute ratio
s_pre = rng;              % pseudorandomly precoder
s_com = rng;              % pseudirabdinly combiner
N_Q = 2;                  % quantization bits
% make dictionary matrices
% load AR_dic and AT_dic Psi_dic from dictionary.mat
% load('dictionary.mat');
AR_dic = Array_dictionary(N_r, G_r);
AT_dic = Array_dictionary(N_t, G_t);
Psi_dic = kron(conj(AT_dic),AR_dic);

%initialize
phi_l   = zeros(sn_length,L);
theta_l = zeros(sn_length,L);
gain    = zeros(L,K,sn_length);
NMSE_dB_times = zeros(sn_length, times);


parfor t = 1:times
    for sn = 1 : sn_length
        %% generate transmit pilot symbols assumed QPSK system
        % symbol_tr(Lt,1,M,K)
%         symbol_tr = generate_symbol(L_t, M, K);
        symbol_tr = pskmod(unidrnd(4,K,1)-1, 4, pi/4, 'gray');
        q_tr = pskmod(unidrnd(4,N_s,M)-1, 4, pi/4, 'gray');

        %% generate channel
        [H_freq, PL] = generate_losnlos_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L, pdr);
%         [H_freq, phi_l, theta_l, alpha, PL] = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);
        % genie aided generate channel
%         [H_freq, phi_l, theta_l, alpha, PL, T_set_genie] = genie_generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L, G_r, G_t);

        %% generate precoder and combiner
        % make training vectors set, A_set = (1,N_Q)
        % make precoder F_tr(Nt,Lt,m) and combiner W_tr(Nr,Lr,m) matrix of m-th training flame
        A_set = 2*pi*(0:N_Q^2-1)/2^(N_Q);
        F_tr = generate_precoder(N_t, L_t, M, A_set, N_Q, s_pre);
        W_tr = generate_combiner(N_r, L_r, M, A_set, N_Q, s_com);

        %% generate measurement matrix Phi
        q_tr_3D = permute(q_tr, [3 1 2]);
        qF_tr = Multiplication3D(q_tr_3D,permute(F_tr, [2 1 3]));
        Phi = zeros(M*L_r, N_t*N_r);
        for m = 1 : M
            Phi((m-1)*L_r+1:m*L_r, :) = kron(qF_tr(:, :, m),W_tr(:, :, m)');
        end

        %% generate covariance matrix and sqrt of covariance matrix
        C_w_m = Multiplication3D(permute(conj(W_tr), [2 1 3]), W_tr);
        C_w = blockdiag(C_w_m);
        D_w = sqrtm(C_w);

        %% caluculate vectorized received signal
        y_k = Phi * reshape(H_freq, N_t*N_r, K) + D_w*crandn(L_r*M, K);
        y = reshape(y_k, L_r*M*K, 1);

        %% calc equivalent measure matrix Upsilon(ML_r*G_tG_r)
        Upsilon = Phi*Psi_dic;

        %% caluculate inverse matrix of D_w
        % D_w_inv can be seen as a frequency flat baseband combiner W_{BB,tr}^(m) 
        % used in the m-th training step
        D_w_inv = inv(D_w);

        %% compute the whitened equivalent obserbation matrix
        % calc whitened equivalent measure matrix Upsilon_w(ML_r*G_tG_r)
        Upsilon_w =(inv(D_w))'*Upsilon;
        

        %% below here, channel estimation
        %% SW-OMP phase
        % initialize    
        MSE_value = 100;
        eta = 1;
        [xi_tilde, T_set] = SW_OMP(MSE_value, y_k, Upsilon_w, D_w, eta, K, M, L_r);
        % genie aided SW-OMP
%         [xi_tilde, T_set] = genie_SW_OMP(y_k, Upsilon_w, D_w, T_set_genie);

        %% LS pahse
%         a_R = sqrt(1/N_r)*exp(1i*2*pi*(0:N_r-1)'.*cos(phi_l_sn));
%         a_T = sqrt(1/N_t)*exp(1i*2*pi*(0:N_t-1)'.*cos(theta_l_sn));
%         Psi = krp(conj(a_T), a_R);
%         Upsilon = Phi * Psi;
%         Upsilon_w =(inv(D_w))'*Upsilon;
%         y_w = (inv(D_w))'*y_k;
%         xi_tilde = pinv(Upsilon)*y_k;        
%         gain_est = reshape(xi_tilde, size(xi_tilde,1)*K,1);
%         krp_array_est = kron(eye(K),krp(conj(a_T),a_R));
%         vec_h = krp_array_est * gain_est;
%         H_freq_est_2 = reshape(vec_h, N_r, N_t, K);
       

        %% channel estimation phase
        [H_freq_est_2,gain_est] = channel_estimation(xi_tilde, K, N_r, N_t, G_t, G_r, T_set);
        % genie channel estimation
%         [H_freq_est_2,gain_est] = genie_channel_estimation(xi_tilde, K, N_r, N_t, G_t, G_r, T_set);

        %% caluculate NMSE at dB domain according to SNR
        NMSE_dB_times(sn, t) = calc_NMSE(H_freq_est_2, H_freq, K, N_r, N_t)-pow2db(N_r*N_t/PL*10^(3));
%         NMSE_dB_times(sn, t) = calc_NMSE(H_freq_est_2, H_freq, K, N_r, N_t)-mean_PL_factor;

    end
end

%% caluculate NMSE at dB domain
NMSE_dB = (sum(NMSE_dB_times,2)/times).'; 
% NMSE_dB_known_angle = (sum(NMSE_dB_times,2)/times).'; 

%% save data
% NMSE_dB_meanPL = NMSE_dB;
NMSE_dB_losnlos = NMSE_dB;
% save("result_SW_OMP_M80.mat", "sn_dB", "NMSE_dB", "times", "M");
% save("result_SW_OMP_M80_losnlos.mat", "sn_dB", "NMSE_dB_losnlos", "times", "M");
% 
%% graph drawing phase
% load("result_SW_OMP_M80.mat");
% save("result_SW_OMP_LosNLos.mat","NMSE_dB_losnlos", "NMSE_dB", "sn_dB", "M", "times");
% figure
% grid on
% hold on
% 
% plot(sn_dB, NMSE_dB, "b-o", "LineWidth", 2, "MarkerSize", 10)
% plot(sn_dB, NMSE_dB_losnlos, "r-s", "LineWidth", 2, "MarkerSize", 10)
% 
% xlabel("SNR (dB)")
% ylabel("NMSE (dB)")
% legend(["SW-OMP (only NLoS)", "SW-OMP (LoS & NLoS)"], "Location","northeast")
% % title(['NMSE versus SNR (M = ', num2str(M), ')'])