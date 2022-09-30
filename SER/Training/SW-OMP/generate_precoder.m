function F_tr = generate_precoder(N_t, L_t, M, A_set, N_Q, s_pre)

    % generate precoder F_tr (N_t, L_t, M)
    rng(s_pre);
    % f_tr_mvec vectorize f_tr for M training step
%     f_tr_mvec = A_set(randi(2^(N_Q),1,N_t*L_t*M))/N_t;
    f_tr_mvec = A_set(randi(2^(N_Q),1,N_t*L_t*M));
    F_tr = reshape(f_tr_mvec,N_t,L_t,M);
    
end