%Generate Problems
function [At,b,c,X,A,C]= Generate_SDP_Problems(n,m)
    %n: n by n symmetric matrices
    %m: number of linear constraints
    rng(1);
    Struc= ones(n);
    %fixed X
    X = 10*sprandsym(Struc);
    eigval = eig(X);
    if (min(eigval)<0)%%make X feasible
        %test = (min(eigval)-1)*speye(size(X));
        X = X - (min(eigval)-1)*speye(size(X));
    end
    eig(X)
    num= n*(n+1)/2;%number of variables
    At = zeros(m,num);
    A = cell(m,1);
    for j=1:m
        A{j} = sprandsym(Struc);
        At(j,:)= ConvMat2LinerCons(A{j});
    end
    Skew_X = X + tril(X,-1);
    %X = X + tril(X,-1);
    idx = Subscribpts(n);
    x = Skew_X(idx)';
    b = At*x;
    
    y = rand(m,1);
    S = 10*sprandsym(Struc);
    eigval = eig(S);
    if (min(eigval)<0)%%make S PSD
        test = (min(eigval)-1)*speye(size(S));
        S = S - (min(eigval)-1)*speye(size(S));
    end
    
    C = S;
    for j = 1:m
        C = C + A{j}*y(j);
    end
    
    %make c a vector
    temp = C+tril(C,-1);
    c = temp(idx)';
end