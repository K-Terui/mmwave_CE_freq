function MSE_dB =  caluculate_MSE_dB(true_value, est_value)
    % caluculate mean square error(MSE)
    % input  : true_value(true value), est_value(estimate value)
    % output : MSE (value of MSE) 
    

    true_vec = reshape(true_value, [], 1);
    est_vec  = reshape(est_value, [], 1);
    N = size(true_vec, 1);
    
    temp = (est_vec - true_vec).*(est_vec - true_vec);
    MSE_dB = pow2db(sum(temp)/N);
    
end