function variance = caluculate_combined_noise_variance(noise_c, M, K, W_tr)

    variance=reshape(var(noise_c),M,K);
    
    for m = 1:M
        theory(:,:,m) = 1*W_tr(:,:,m)'*W_tr(:,:,m);
    end
    
end