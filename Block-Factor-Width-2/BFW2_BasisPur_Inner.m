%Basis Pursuit Block Factor Width 2
%Inner Approximation

function [x,y,info,OBJ] = BFW2_BasisPur_Inner(A,b,C,P)
%Solve in Sedumi "Primal" form
%Input:
%A: a cell storing a sequence of symmetric matrices.
%b: a column vector.
%C: the coefficient matrix of objective function    
    
    dim_mat = width(C);%dimension of design matirx
    m = length(A);%number of linear constraints

    %check partition 
    if (sum(P)~=dim_mat)
        fprintf("The size of partition does not match the dimension of decision matrix");
        return ; 
    end
    
    %construct Operator (substract sub-principle matrix)
    % combination: 
    CardOfP = length(P);
    Comb = nchoosek(1:CardOfP,2);%block factor-width-2
    I = eye(dim_mat);
    E = cell(1,CardOfP);
    
    %define a set of integers (M)
    CardOfM = CardOfP+1;
    M = zeros(1,CardOfM);
    M(1) = 1;
    for i = 2:CardOfM
        M(i) = M(i-1)+P(i-1);
    end
    count = 0; %count number of variables
    for i = 1:CardOfP
        E{i} = I(M(i):(M(i+1)-1),:);
        count = count + width(E{i})^2;
    end
    
    U = eye(dim_mat);

    %In each iteration, Objection function and constraints change
    CNVG = false;
    %Convergence condition happens when abs(PreObj-NowObj)) < err
    Iter = 1; err = 1.0e-3;
    OBJ = []; %keep track of Obj Value
    PreX = zeros(dim_mat); 
    %epsilon = 0.00001;
    
    while ~CNVG
                
        %form Obj function
        c_new = [];
        A_new = [];
        K.s = [];

        NumOfComb = length(Comb);
        if length(P) <=2
            NumOfComb = 1;
        end

        for n = 1:NumOfComb
            i = Comb(n,1); j = Comb(n,2);
            EE = U'*([E{i};E{j}]');
            c_new = [c_new;vec(EE'*C*EE)];
            [numRows,numCols] = size(EE);
            %if (Iter == 1)
                K.s = [K.s numCols];
            %end
        end
            
        %constraint
        for k = 1:length(A)
            temp = [];
            for n=1:NumOfComb
                i = Comb(n,1); j = Comb(n,2);
                EE = U'*([E{i};E{j}]'); 
                [numRows,numCols] = size(EE);
                temp = [temp;vec(EE'*A{k}*EE)];
            end
            A_new = [A_new, temp];
        end   
        
        [x,y,info] = sedumi(A_new',b,c_new,K);
    
        %Recover X from x
        X = zeros(length(C));
        s = 1;
        for n = 1:NumOfComb
             i = Comb(n,1); j = Comb(n,2);
             EE = U'*([E{i};E{j}]');
             tempx = x(s:s+K.s(n)^2-1);
             tempx = reshape(tempx,[K.s(n) K.s(n)]);
             s = s+K.s(n)^2;
             X = X + EE*tempx*EE';
        end
        
        OBJ = [OBJ, trace(C'*X)];
        
        if (Iter>1)
            if  (OBJ(Iter)-OBJ(Iter-1))^2<err
                CNVG = true;
                break;
            end
        end
        PreX = X;
        Iter = Iter+1;
        D = eig(PreX)
        U = chol(PreX);
        
    end
end