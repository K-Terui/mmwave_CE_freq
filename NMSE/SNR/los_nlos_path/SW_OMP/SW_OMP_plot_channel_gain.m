%{
SW-OMP NMSE
%}
% make received signal
clear
% definition
times = 1280;             % Monte Carlo times
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
% x_tilde_store = zeros(G_r*G_t,1,K);
% MSE_store = zeros(sn_length);
NMSE_dB = zeros(1, sn_length);
store_fro_est = zeros(sn_length, K, times);
est_channel_gain_store = zeros(sn_length, K, times);

parfor t = 1:times
    for sn = 1 : sn_length
    %%
    % make transmit pilot assume QPSK system
    % symbol_tr(Lt,1,M,K)
    symbol_tr = generate_symbol(L_t, M, K);

    %%
    % generate channel
    % initialize matrices
%     a_R=zeros(N_r,1);
%     a_T=zeros(N_t,1);
%     R=zeros(1,N_r);
%     T=zeros(1,N_t);
    noise_c=zeros(L_r,1,M,K);
    y_received=zeros(L_r,1,M,K);

    % make H_freq[k] (Nr,Nt,k) (Heath eq(2)) and noise^(m)[k] (Nr,1,m,k)
    [H_freq, phi_l, theta_l, gain] = generate_channel(K, N_C, N_r, N_t, sn_dB(sn), power_t, L);
    noise = generate_noise(K, M, N_r);

    %%
    % make precoder and combiner
    % make precoder F_tr(Nt,Lt,m) and combiner W_tr(Nr,Lr,m) matrix of m-th training flame
    % make training vectors set, A_set = (1,N_Q)
    A_set = 2*pi*(0:N_Q^2-1)/2^(N_Q);
    F_tr = generate_precoder(N_t, L_t, M, A_set, N_Q, s_pre);
    W_tr = generate_combiner(N_r, L_r, M, A_set, N_Q, s_com);
    
    %%
    % noise combining n_c^(m)[k] noise_c(Lr,1,m,k)
    for k = 0:K-1
        for m = 1:M
            noise_c(:,:,m,k+1) = W_tr(:,:,m)'*noise(:,:,m,k+1); 
        end
    end

    %%
    % make received signal y[k] y_received(Lr,1,M,K)
    for k = 0:K-1
        for m = 1:M
            y_received(:,:,m,k+1) = ctranspose(W_tr(:,:,m))*H_freq(:,:,k+1)*F_tr(:,:,m)*symbol_tr(:,:,m,k+1)+noise_c(:,:,m,k+1);
        end
    end

    %%
    % make vectorized measurement matrix Phi_measure(MLr,NtNr)
    % initialize matrices
    Phi_measure_mth = zeros(L_r, N_t*N_r,M);
    Phi_measure = zeros(M*L_r, N_t*N_r);
    for m = 1:M
        Phi_measure_mth(:,:,m) = 1/sqrt(2)*kron(transpose(F_tr(:,:,m)),ctranspose(W_tr(:,:,m)));
        if m == 1
            Phi_measure = Phi_measure_mth(:,:,m);
        else
            Phi_measure = vertcat(Phi_measure,Phi_measure_mth(:,:,m));
        end
    end

    %%
    % make vectorized combined noise noise_vec(MLr,1,k) from noise_combiner(Lr,1,m,k)
    % make vectorized received signal y_received_vec(MLr,1,k)
    noise_vec = zeros(M*L_r,1,K);
    noise_vec_kth = zeros(M*L_r,1);
    y_received_vec = zeros(M*L_r,1,K);
    y_received_vec_kth = zeros(M*L_r,1);
    for k = 0:K-1
        for m = 1:M
            if m == 1
                noise_vec_kth = noise_c(:,:,m,k+1);
                y_received_vec_kth = y_received(:,:,m,k+1);
            else
                noise_vec_kth = vertcat(noise_vec_kth,noise_c(:,:,m,k+1));
                y_received_vec_kth = vertcat(y_received_vec_kth,y_received(:,:,m,k+1));
            end
        end
        noise_vec(:,:,k+1) = noise_vec_kth;
        y_received_vec(:,:,k+1) = y_received_vec_kth;
    end

    %%
    % calc equivalent measure matrix Upsilon(ML_r*G_tG_r)
    Upsilon = Phi_measure*Psi_dic;

    %%
    % noise covariance matrix of y[k] is C_w
    % calc C_w=blkdiag{W^(1)^H*W^(1),...,} C_w(MLr,MLr)
    for m = 1:M
        if m == 1
            C_w = W_tr(:,:,m)'*W_tr(:,:,m);
        else
            C_w = blkdiag(C_w,W_tr(:,:,m)'*W_tr(:,:,m));
        end
    end

    %%
    % Cholesky factorization  C_w = D_w^H D_w
    % D_w(MLr,MLr) this is upper triangular matrix
    D_w = chol(C_w);

    %%
    % calc D_w_inv(MLr,MLr)
    % D_w_inv can be seen as a frequency flat baseband combiner W_{BB,tr}^(m) 
    % used in the m-th training step
    D_w_inv = inv(D_w);

    %%
    % calc whitened equivalent measure matrix Upsilon_w(ML_r*G_tG_r)
    Upsilon_w = ctranspose(D_w_inv)*Upsilon;

    %%
    % calc  correlation vector cor_vec (GtGr,1,K)
    % cor_vec = (D_w_inv'*Upsilon)'*D_w_inv'*y_received_vec;
    cor_vec = zeros(G_t*G_r,1,K);
    for k = 0:K-1
        cor_vec(:,:,k+1) = Upsilon_w'*D_w_inv'*y_received_vec(:,:,k+1);
    end

    % %%
    % % compute the whitende equivalent observation matrix
    % Upsilon_w = D_w_inv'*Phi_measure*Psi_dic;

    %% SW-OMP
    MSE_value = 100;
    eta = 1;
    [xi_tilde, T_set] = SW_OMP(MSE_value, y_received_vec, Upsilon, D_w, eta, G_r, G_t, K, M, L_r);

    %% Channel Estimation phase
    [H_freq_est_2, phi_est_sn, theta_est_sn, gain_est_sn] = channel_estimation_2(xi_tilde, T_set, K, N_r, N_t, AR_dic, AT_dic, G_t, G_r);
    
    %% caluculate estimated channel gain
    est_vec = reshape(H_freq_est_2, N_r*N_t, K);
    fro_est = vecnorm(est_vec).^2;
    store_fro_est(sn,:,t) = pow2db(fro_est);

    end
end

est_channel_gain = reshape(store_fro_est, sn_length, K*times);

%% plot a glaph 
% 
% hold off
% figure
% hold on
% 
% for sn = 1 : sn_length 
%     
%     cdf_fig = cdfplot(est_channel_gain(sn,:));
%     cdf_fig.LineWidth = 2;
%     
% end
% 
% xlabel("Channel Gain (dB)")
% ylabel("Comulative Distribution")
% legend(["SNR = -15dB", "SNR = -10dB","SNR = -5dB", "SNR = 0dB", "SNR = 5dB", "SNR = 10dB",], "Location","southeast")
% % title(['CDF (K = ', num2str(K), ')'])
% xlabel("Channel Gain (dB)")
% ylabel("Comulative Distribution")
% legend(["SNR = -15dB", "SNR = -10dB","SNR = -5dB", "SNR = 0dB", "SNR = 5dB", "SNR = 10dB",], "Location","southeast")
% % title(['CDF (K = ', num2str(K), ')'])

