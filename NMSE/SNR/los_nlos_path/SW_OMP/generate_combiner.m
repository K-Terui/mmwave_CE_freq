function W_tr = generate_combiner(N_r, L_r, M, A_set, N_Q, s_com)

    % generate combiner W_tr (N_r, L_r, M)
    rng(s_com);
    % w_tr_mvec vectorize w_tr for M training step
%     w_tr_mvec = A_set(randi(2^(N_Q),1,N_r*L_r*M))/N_r;
    w_tr_mvec = A_set(randi(2^(N_Q),1,N_r*L_r*M));
    W_tr = reshape(w_tr_mvec,N_r,L_r,M);
    
end