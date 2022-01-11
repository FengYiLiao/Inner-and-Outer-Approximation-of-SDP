%Basis Pursuit Block Factor Width 2
%Inner Approximation

function [x,y,info,OBJ] = BFW2_BasisPur_Inner(A,b,C,P)
%Solve in Sedumi "Primal" form
%Input:
%A: a cell storing a sequence of symmetric matrices.
%b: a column vector.
%C: the coefficient matrix of objective function    
    Max_Iter = 200;
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
    NumOfComb = length(Comb);
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
    
%     %Test
%         %Outer Approximation
%         %Solve in Sedumi dual form
%         b_Out = -vec(C);%minimize
%         c_Out = b; %linear constraints %will add more later on
%         K_Out.f = m;
%         K_Out.s=[];
% 
%         %Outer-Linear constraints
%         At_Out = [];
%         for i = 1:m% m is number of constraints
%             At_Out = [At_Out; vec(A{i})'];
%         end
% 
%         %Outer-symmetric constraints
%         temp = ones(dim_mat);
%         temp = tril(temp,-1);
%         idx = find(temp); %index of symmetric element
%         for i = 1:length(idx)
%             temp_mat = zeros(dim_mat);
%             temp_mat(idx(i)) = 1;
%             temp_mat = temp_mat + (-temp_mat');
%             At_Out = [At_Out; vec(temp_mat)'];
%         end
%         c_Out = [c_Out;zeros(length(idx),1)];%Symmetric Constraints
% 
%         K_Out.f = m+length(idx);%Linear Constraints and Symmetric Constraints
% 
%         %Outer-Cone constraints
%         %In each iteration, the only one change happens here
%         TempAt=[];
%         count = 0; %count number of variables
%         for num = 1:NumOfComb
%             i = Comb(num,1); j = Comb(num,2);
%             EE = [E{i};E{j}]*U;
%             tempAt = [];
%             for r1 = 1: width(EE)
%                 for r2 =1:width(EE)
%                     tempAt = [tempAt, -vec(EE(:,r1)*EE(:,r2)')];
%                 end
%             end
%             count = count + height(EE)^2;
%             K_Out.s = [K_Out.s height(EE)];
%             TempAt=[TempAt;tempAt];
%         end    
%         At_Out = [At_Out;TempAt];
%         At_Out = At_Out';
%         c_Out = [c_Out;zeros(count,1)];
%     
%         [x_Out,y_Out,info_Out]=sedumi(At_Out,b_Out,c_Out,K_Out);
%         
%         %projection
%         X_Out = mat(y_Out); X_Out = full(X_Out);
%         [EigVec,EigVal] = eig(X_Out);
%         EigVal(EigVal<0) = 0;
%         X_Out = EigVec*EigVal*EigVec';
%         
%         [V D]=eig(X_Out); V=real(V);D=real(D);
%         Norm = sqrt(sum(V.*V));
%         V = V./(ones(dim_mat,1)*Norm);
%         U = (V*sqrt(D))'; U=real(U);
%     %Test
    
    
    

    %In each iteration, Objection function and constraints change
    CNVG = false;
    %Convergence condition happens when abs(PreObj-NowObj)) < err
    Iter = 1; err = 1.0e-3; Count = 0; 
    OBJ = []; %keep track of Obj Value
    PreX = zeros(dim_mat); 
    %epsilon = 0.00001;
    epsilon = 0.001;
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
            if  (OBJ(Iter)-OBJ(Iter-1))^2<err || Iter>=Max_Iter
                Count = Count +1;
                if (Count >= 3)
                    CNVG = true;
                    break;
                end
            else 
                Count = 0;
            end
            
        end
        PreX = X;
        Iter = Iter+1;
        D = eig(PreX);
        if (sum(D<0)>0)%address numerical issue
            PreX = PreX + abs(min(D)+0.00001)*eye(dim_mat);
        end
        
%Change here for different decomposition
%Cholesky Decomposition
%         D = eig(PreX);
%         if (min(D)<=10^(-10))
%             PreX = PreX + 10^(-10)*eye(dim_mat);
%         end
%        D = eig(PreX)
%        U = chol(PreX);

%try different decompostion, which works for PSD matrices
%Square Root Decomposition
        [V D]=eig(PreX);
        Norm = sqrt(sum(V.*V));
        V = V./(ones(dim_mat,1)*Norm);
        U = (V*sqrt(D))';
        
        
        
    end
end
