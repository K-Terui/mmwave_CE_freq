function Y = pagectranspose(X)
    
    % page wise conjugate transpose matrix
    Y = permute(conj(X),[2 1 3:ndims(X)]);
    
end