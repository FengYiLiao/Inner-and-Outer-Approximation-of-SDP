function g = ConvMat2LinerCons(B)
    %Input the new atom B
    %Output the linear constraint
    [rows cols] = size(B);
    %NewB = B + tril(B,-1);%because of symmetric and linear constraint
    NewB = B+triu(B,1)'; %Add Upper triangle
    n = rows*(rows+1)/2;%dimension
    g = zeros(n,1);
    LB= tril(NewB);%lower triangle of B
    [row col]=find(LB);%Find indices and values of nonzero elements
    index = (row-1).*row/2+col; %index of the matlab matrix format
    for i = 1 : length(row)
        g(index(i)) = NewB(row(i),col(i));
    end
end
