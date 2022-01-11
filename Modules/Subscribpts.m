function ind = Subscribpts(k)
    %k x k symmetric matrix
    %get the indices of the decision elements
    row =[];
    col =[];
    for j = 1 : k
        row = [row [j*ones(1,j)]];
        col = [col [1:j]];
    end
    ind = sub2ind([k k],row,col);
end
