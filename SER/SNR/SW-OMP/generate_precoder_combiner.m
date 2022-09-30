function [F_tr, W_tr] = generate_precoder_combiner(N_t, N_r, L_t, L_r, M, A_set, N_Q, s_pre, s_com, H_freq, K)
%     %% generate precoder F_tr (N_t, L_t, M)
%     rng(s_pre);
%     % f_tr_mvec vectorize f_tr for M training step
% %     f_tr_mvec = A_set(randi(2^(N_Q),1,N_t*L_t*M))/N_t;
%     f_tr_mvec = A_set(randi(2^(N_Q),1,N_t*L_t*M));
%     F_tr = reshape(f_tr_mvec,N_t,L_t,M);
%     
%     %% generate combiner W_tr (N_r, L_r, M)
%     rng(s_com);
%     % w_tr_mvec vectorize w_tr for M training step
% %     w_tr_mvec = A_set(randi(2^(N_Q),1,N_r*L_r*M))/N_r;
%     w_tr_mvec = A_set(randi(2^(N_Q),1,N_r*L_r*M));
%     W_tr = reshape(w_tr_mvec,N_r,L_r,M);

%% Precoder and Combiner design from SVD
    for k = 0:K-1
        [W_tr(:,:,k+1), S, F_tr(:,:,k+1)] = svd(H_freq(:,:,k+1));
    end

end