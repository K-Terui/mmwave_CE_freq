%{
OMP NMSE
%}
% make received signal
clear
% definition
times = 1000;             % Monte Carlo times
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
% make dictionary matrices
% load AR_dic and AT_dic Psi_dic from dictionary.mat
% load('dictionary.mat');
AR_dic = Array_dictionary(N_r, G_r);
AT_dic = Array_dictionary(N_t, G_t);
Psi_dic = kron(conj(AT_dic),AR_dic);

%initialize
% x_tilde_store = zeros(G_r*G_t,1,K);
% MSE_store = zeros(length(sn_dB));
NMSE_dB = zeros(1, length(sn_dB));


for sn = 1 : length(sn_dB)
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
%     noise = zeros(N_r,1,M, K);

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
    % OMP phase
    % initialize
    MSE_value = 100;
    x_hat = zeros(G_r*G_t, K);
    index = zeros(L,K);

    for k = 0:K-1

        A = Upsilon;
        A_hat = zeros(size(A));
        I = eye(G_r*G_t);
        r = y_received_vec(:, :, k+1);

        % OMP minimize x_hat
        % subject to |y-Ax|^2_2<eta

        for i = 1:L
            % decide index
            [~, index(i, k+1)] = max(abs(A'*r));

            % obtain column vector from measurement matrix
            A_hat(:, i) = A(:, index(i, k+1));

            % clean up matrix A
            A(:, index(i, k+1)) = 0;

            % LS
            x_tilde = A_hat(:, 1:i)\y_received_vec(:, :, k+1);

            % apdate residual
            r = y_received_vec(:, :, k+1) - A_hat(:, 1:i)*x_tilde;
        end
        % reconstruct sparse signal
        x_hat(index(:,k+1), k+1) = x_tilde;
    end
 
%%
    % Channel Estimation phase
    H_freq_est = channel_estimation(x_hat, index, K, N_r, N_t, AR_dic, AT_dic, phi_l, theta_l, alpha_l, G_r, G_t, N_C);
    NMSE_dB(1,sn) = calc_NMSE(H_freq_est, H_freq, K, N_r, N_t);
%     x_hat_store(:,:,:,sn) = index;
%     MSE_value_store(sn) = MSE_value;
end
%%
% % save data
% save("result.mat", "sn_dB", "NMSE_dB");
% 
%% 
% % plot a glaph 
% load("Theory.mat");
% figure
% grid on
% hold on
% plot(sn_dB, NMSE_dB, "b-o", "LineWidth", 2, "MarkerSize", 10)
% plot(sn_dB, NMSE_theory_OMP, "g-s", "LineWidth", 2, "MarkerSize", 10)
% 
% xlabel("SNR (dB)")
% ylabel("NMSE (dB)")
% legend(["OMP (Simu)", "OMP (Conv)"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])
% 
% figure
% plot(sn_dB, NMSE_dB, "b-o", "LineWidth", 2, "MarkerSize", 10)
% grid on
% xlabel("SNR (dB)")
% ylabel("NMSE (dB)")
% legend("OMP", "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])



