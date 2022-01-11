function [x,y,info,OBJ,X] = BFW2(A,b,C,P,U)
%Solve in Sedumi "Primal" form
%Input:
%     A: a cell storing a sequence of symmetric matrices.
%     b: a column vector.
%     C: the coefficient matrix of objective function
%     U: rotation matrix for extreme rays. 
%Output:
%     X: decision matrix 

%Use Factor width-2 to solve the problem
%P: partition. sum(P) = n (for a n X n matrix)
%Note: only does one iteration

    c = vec(C);
    dim_mat = sqrt(length(c));
    
    %construct Operator (substract sub-principle matrix)
    % combination: 
    CardOfP = length(P);
    Comb = nchoosek(1:CardOfP,2);%block factor-width-2
    I = eye(dim_mat);
    E = cell(1,CardOfP);
    
    %define a set of integers (m)
    CardOfm = CardOfP+1;
    m = zeros(1,CardOfm);
    m(1) = 1;
    for i = 2:CardOfm
        m(i) = m(i-1)+P(i-1);
    end
    for i = 1:CardOfP
        E{i} = I(m(i):(m(i+1)-1),:);
    end
    
    %form Obj function
    c_new = [];
    At_new = [];
    K.s = [];
    
    NumOfComb = length(Comb);
    if length(P) <=2
        NumOfComb = 1;
    end
    
    for n = 1:NumOfComb
        i = Comb(n,1);
        j = Comb(n,2);
        
        %EE = [ E{i};E{j}];
        %c_new = [c_new;vec(EE*C*EE')];
        
        EE = U'*([E{i};E{j}]');
        c_new = [c_new;vec(EE'*C*EE)];
        
        [numRows,numCols] = size(EE);
        K.s = [K.s numCols];
    end
    
    %constraint
    for k = 1:length(A)
        temp = [];
        for n=1:NumOfComb
            i = Comb(n,1);
            j = Comb(n,2);
            
%             EE = [ E{i};E{j}]; 
%             [numRows,numCols] = size(EE);
%             temp = [temp;vec(EE*A{k}*EE')];
            
            
            EE = U'*([E{i};E{j}]'); 
            %[numRows,numCols] = size(EE);
            temp = [temp;vec(EE'*A{k}*EE)];
        end
        At_new = [At_new, temp];
    end 
    
    [x,y,info] = sedumi(At_new',b,c_new,K);
    
    %Recover X from x
    X = zeros(length(C));
    s = 1;
    for n = 1:NumOfComb
         i = Comb(n,1);
         j = Comb(n,2);
         
%          EE = [E{i};E{j}];
%          tempx = x(s:s+K.s(n)^2-1);
%          tempx = reshape(tempx,[K.s(n) K.s(n)]);
%          s = s+K.s(n)^2;
%          X = X + EE'*tempx*EE;
         
         
         EE = U'*([E{i};E{j}]');
         tempx = x(s:s+K.s(n)^2-1);
         tempx = reshape(tempx,[K.s(n) K.s(n)]);
         s = s+K.s(n)^2;
         X = X + EE*tempx*EE';
    end
    OBJ = trace(C'*X);

end
