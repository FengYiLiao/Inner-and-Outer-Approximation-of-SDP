%LP column generation using Sedumi
%Sedumi dual standard form
%Outer Approximation
%This is the first version of LPColsGen. The format is less
%straightforward. The format is different from Basis Pursuit
%Index of elements in the design matrix follows below
%X is n x n
%X = [x1 ..................
%     x2 x3 ............... 
%     x4 x5 x6  ........... 
%     x7 x8 x9 x10.........
%         .       . 
%         .         .
%         .           .
%     .................x(n*(n+1)/2)]  
function [x,y,info,OBJ]=LPColsGen(At,b,c,K)
%input:
%       At, b, c, K. 
%       At is a matrix. b is a vector. c is a vector.
%       solve in Sedumi dual form   
%Output:
%       x,y,info,OBJ 
%       y constains the decision variables we need
%       OBJ consains the history of the Obj value
    MaxItr = 200;
    At_new = -At';
    [M N]=size(At_new);
    dim_mat = (-1+sqrt(1+8*M))/2; %dimension of the symmetric matrix
    idx = Subscribpts(dim_mat);%get the indices of the desgin elements in the matrix
    b_new = -c;
    c_new = -b;
    K_new = K;
    %check the cone
    if ~isfield(K_new,'l')
        K_new.l = 0;
    end
    if ~isfield(K_new,'q')
        K_new.q = [];
    end
    
    %%add fixed atoms %rank one matrix
    Comb = nchoosek(1:dim_mat,2);
    NumOfComb = length(Comb);
    if dim_mat <3
        NumOfComb = 1;
    end
    Atom = zeros(dim_mat);%new Atom
    for j = 1:NumOfComb
        for i = 1:2
            u = zeros(1,dim_mat)';
            if (i == 1)
                u(Comb(j,1))=1;
                u(Comb(j,2))=-1;
            else
                u(Comb(j,1))=1;
                u(Comb(j,2))=1;
            end
            Atom = u*u';
            g = ConvMat2LinerCons(Atom); %Convert matrix into linear constraint following the indexing format 
            At_new = [At_new,-g];
            c_new = [c_new; 0];
            K_new.l= K_new.l+ 1;
        end
    end

    
    PSD = false;
    Iter = 1;
    err = 1.0e-3;
    OBJ = [];
    while ~PSD
        % Compute the optimization.
        [x,y,info]=sedumi(At_new,b_new,c_new,K_new);
        OBJ = [OBJ, -b_new'*y]; 
        X = zeros(dim_mat);
        X(idx) = y;% Show the optimal x solution.
        X = X + tril(X,-1)';
        [V,D] = eig(X);
        [d, I]=sort(diag(D));
        if d(1)>=0 || abs(d(1))<=err || Iter>MaxItr %get PSD 
            PSD = true;
            break;
        end
        u = V(:,I(1));%the most negative eigenvector
        Atom = u*u';%new Atom
        g = ConvMat2LinerCons(Atom);%linear constraint
        At_new = [At_new,-g];%one more linear constraint
        c_new = [c_new; 0];%one more linear constraint
        K_new.l= K_new.l+ 1;%one more linear constraint
        Iter = Iter+1;
    end
    disp("Iter")
    disp(Iter);
end
