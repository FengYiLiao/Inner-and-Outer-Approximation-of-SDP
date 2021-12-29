%SOCP column generation using Sedumi
%dual standard form
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
function [x,y,info,OBJ]= SOCPColsGen(At,b,c,K)
%input:
%       At, b, c, K. 
%       At is a matrix. b is a vector. c is a vector.
%       solve in Sedumi dual form   
%Output:
%       x,y,info,OBJ 
%       y constains the decision variables we need
%       OBJ consains the history of the Obj value
    MaxItr = 200;
    At_new = sparse(-At');
    [M N]=size(At_new);
    dim_mat = (-1+sqrt(1+8*M))/2; %dimension of the symmetric matrix
    idx = Subscribpts(dim_mat);%get the indices of the desgin elements in the matrix
    b_new = sparse(-c);
    c_new = sparse(-b);
    K_new = K;
    %check the cone
    if ~isfield(K_new,'l')
        K_new.l = 0;
    end
    if ~isfield(K_new,'q')
        K_new.q = [];
    end
 
    dim_Var = length(b_new);
    
    %add fixed atoms
    
    Comb = nchoosek(1:dim_mat,2);
    NumOfComb = length(Comb);
    if dim_mat <3
        NumOfComb = 1;
    end
    At_l = sparse(dim_Var,2*NumOfComb);%linear constraints% zeros(dim_Var,2*NumOfComb);
    At_q = sparse(dim_Var,3*NumOfComb);%quadratic constraints% zeros(dim_Var,3*NumOfComb);
    l= 1; q = 1;
    for j = 1:NumOfComb
        V = zeros(dim_mat,2);
        V(Comb(j,1),1)=1;
        V(Comb(j,2),2)=1;
        [C1 C2 C3] = ConvMat2SOCPCons(V,dim_Var,dim_mat);%Convert matrix into linear constraint following the indexing format 
        
        At_l(:,l) = -C1';
        At_l(:,l+1)=-C3';
        l = l+2;
        At_q(:,q) = -(C1+C3)';
        At_q(:,q+1)=-2*C2';
        At_q(:,q+2)=-(C1-C3)';
        q = q + 3;        
    end
    %update the constraint
    At_new = [At_new(:,1:K_new.f+K_new.l),At_l,At_new(:,K_new.f+K_new.l+1:end),At_q];
    c_new = [c_new;sparse(5*NumOfComb,1)];%zeros(5*NumOfComb,1)
    K_new.l = K_new.l+2*NumOfComb;
    K_new.q = [K_new.q 3*ones(1,NumOfComb)];
       
    PSD = false; Iter = 1; err = 1.0e-3;
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
        v = V(:,I(1:2));
        v(:,1)=v(:,1)/norm(v(:,1));
        v(:,2)=v(:,2)/norm(v(:,2));
        Iter = Iter+1;
        [C1 C2 C3] = ConvMat2SOCPCons(v,dim_Var,dim_mat);%coefiecients of linear constraint in the symmetric matrix
        At_new = [At_new(:,1:K_new.f+K_new.l),sparse(-C1'),sparse(-C3'),At_new(:,K_new.f+K_new.l+1:end),sparse(-(C1+C3)'),sparse(-(C1-C3)'),sparse(-2*C2')];
        c_new = [c_new;sparse(5,1)];%2 linear constraints and 3 quadradic constraints
        K_new.l= K_new.l+ 2;%2 linear constraints
        K_new.q =[K_new.q 3];%3 quadradic constraints
    end
    disp("Iter")
    disp(Iter);
end