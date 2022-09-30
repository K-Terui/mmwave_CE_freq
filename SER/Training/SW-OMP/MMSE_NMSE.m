%{
SW-OMP NMSE
%}
% make received signal
clear
% definition
times = 1;             % Monte Carlo times (9/1 5600)
M = 100;                  % training flame length 
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
% phi_est     = zeros(sn_length,L);
% theta_est   = zeros(sn_length,L);
% gain_est    = zeros(L,K,sn_length);
% x_tilde_store = zeros(G_r*G_t,1,K);
% MSE_store = zeros(sn_length);
NMSE_dB_times = zeros(sn_length, times);

for t = 1:times
    for sn = 1 : sn_length
        %%
        % make transmit pilot assume QPSK system
        % symbol_tr(Lt,1,M,K)
%         symbol_tr = generate_symbol(L_t, M, K);
        symbol_tr = pskmod(unidrnd(4,K,1)-1, 4, pi/4, 'gray');
        q_tr = pskmod(unidrnd(4,N_s,M)-1, 4, pi/4, 'gray');

        %%
        % generate channel
        % initialize matrices
    %     a_R=zeros(N_r,1);
    %     a_T=zeros(N_t,1);
    %     R=zeros(1,N_r);
    %     T=zeros(1,N_t);

        % make H_freq[k] (Nr,Nt,k) (Heath eq(2)) and noise^(m)[k] (Nr,1,m,k)
    %     H_freq = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);
        [H_freq, phi_l_sn, theta_l_sn, gain_sn] = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);
        
%         phi_l(sn,:) =  phi_l_sn;
%         theta_l(sn,:) = theta_l_sn;
%         gain(:,:,sn) = gain_sn;

        %%
        % make precoder and combiner
        % make precoder F_tr(Nt,Lt,m) and combiner W_tr(Nr,Lr,m) matrix of m-th training flame
        % make training vectors set, A_set = (1,N_Q)
        A_set = 2*pi*(0:N_Q^2-1)/2^(N_Q);
        F_tr = generate_precoder(N_t, L_t, M, A_set, N_Q, s_pre);
        W_tr = generate_combiner(N_r, L_r, M, A_set, N_Q, s_com);

        %%
        q_tr_3D = permute(q_tr,[3 1 2]);
        qF_tr = Multiplication3D(q_tr_3D,permute(F_tr,[2 1 3]));
        Phi = zeros(M*L_r,N_t*N_r);
        for m = 1 : M
            Phi((m-1)*L_r+1:m*L_r,:) = kron(qF_tr(:,:,m),W_tr(:,:,m)');
        end
        C_w_m = Multiplication3D(permute(conj(W_tr),[2 1 3]),W_tr);
        C_w = blockdiag(C_w_m);
        D_w = sqrtm(C_w);

        %%
        y_k = Phi * reshape(H_freq,N_t*N_r,K) + D_w*crandn(L_r*M,K);
        y = reshape(y_k,L_r*M*K,1);

        %%
        % calc equivalent measure matrix Upsilon(ML_r*G_tG_r)
        Upsilon = Phi*Psi_dic;

        %%
        % calc D_w_inv(MLr,MLr)
        % D_w_inv can be seen as a frequency flat baseband combiner W_{BB,tr}^(m) 
        % used in the m-th training step
        D_w_inv = inv(D_w);

        %% Compute the whitened equivalent obserbation matrix
        % calc whitened equivalent measure matrix Upsilon_w(ML_r*G_tG_r)
        Upsilon_w =(inv(D_w))'*Upsilon;
        

        %%
        % MMSE phase
        % initialize    
        MSE_value = 100;
        eta = 1;
        [xi_tilde, T_set] = MMSE(MSE_value, Upsilon_w, y_k, D_w, K, M, L_r, eta);
       

    %%
%         % Channel Estimation phase
        [H_freq_est_2, gain_est_sn] = channel_estimation_2(xi_tilde, K, N_r, N_t, AR_dic, AT_dic, G_t, G_r, T_set);
        NMSE_dB_times(sn, t) = calc_NMSE(H_freq_est_2, H_freq, K, N_r, N_t);
%         [H_freq_est_3,gain_est] = channel_estimation_3(xi_tilde, K, N_r, N_t, phi_l_sn, theta_l_sn, G_t, G_r, T_set);
%         NMSE_dB_times(sn, t) = calc_NMSE(H_freq_est_3, H_freq, K, N_r, N_t);
%         

    end
end

NMSE_dB = (sum(NMSE_dB_times,2)/times).';
% %%
% % save data
% save("result.mat", "sn_dB", "NMSE_dB");
% 
%% plot a glaph 

% figure
% plot(sn_dB, NMSE_dB, "b-o", "LineWidth", 2, "MarkerSize", 10)
% xlabel("SNR (dB)")
% ylabel("NMSE (dB)")
% legend("MMSE", "Location","northeast")
% % title(['NMSE versus SNR (M = ', num2str(M), ')'])
% grid on
