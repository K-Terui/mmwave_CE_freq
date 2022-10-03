function noise = generate_noise(K, M, N_r)
    % calc circularly symmetric complex Gausssian 
    % distributed additive noise vector (frequency domain), 'noise'
    % size(noise) = N_r, 1, M, K
    % n[k] ~ CN(0,1)

    vec_n = crandn((N_r*M*K),1);
    noise = reshape(vec_n,N_r,1,M,K);

end