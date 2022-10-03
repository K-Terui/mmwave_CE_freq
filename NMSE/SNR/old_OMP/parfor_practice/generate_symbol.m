function symbol_tr = generate_symbol(L_t, M, K)
    % make transmit symbol symbol_tr (L_t,1,M,K)
    
    x_rand = randi(2,1,M*K);
    y_rand = randi(2,1,M*K);
    x = (x_rand>1).*(-1) + (x_rand==1).*1;
    y = (y_rand>1).*(-1) + (y_rand==1).*1;
    
    vec_symbol = (x+y*1i)/sqrt(2);
    symbol_tr = reshape(vec_symbol,L_t,1,M,K);
    
end