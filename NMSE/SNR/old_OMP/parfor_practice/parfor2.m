parfor x = 0:10
    for y = 2:10
        A(y) = A(y-1) + y;
    end
end