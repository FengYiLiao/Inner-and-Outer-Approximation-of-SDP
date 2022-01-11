%Basis Pursuit Block Factor Width 2 
%Outer Approxima%Basis Pursuit Outer Approximation relies on the inner approximation of
%generating a sequence of rotation matrices, we denote it as U in the
%script. 

function [x_Out,y_Out,info_Out,OBJ_Out,OBJ_U] = BFW2_BasisPur_Outer(A,b,C,P,Mode)
%Solve in Sedumi "Primal" form
%Input:
%A: a cell storing a sequence of symmetric matrices.
%b: a column vector.
%C: the coefficient matrix of objective function    
%Mode: converge condition. 1:PSD  2.Duality gap < epsilon 
    epsilon = 0.001;
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
    I = eye(dim_mat);
    E = cell(1,CardOfP);
    
    %define a set of integers (M)
    CardOfM = CardOfP+1;
    M = zeros(1,CardOfM);
    M(1) = 1;
    for i = 2:CardOfM
        M(i) = M(i-1)+P(i-1);
    end
    
    for i = 1:CardOfP
        E{i} = I(M(i):(M(i+1)-1),:);
    end
    
    U = eye(dim_mat);%first U
    
    %Start with inner approximation
    %Sovle in sedumi primal form

    %Linear Constraints constraints
    A_In = [];
    K_In.f = m;
    K_In.s= [];
    for i = 1:m% m is number of constraints
        A_In = [A_In, vec(A{i})];
    end
    
    %Cone constraints- Right hand side
    %In each iteration, the only one change happens here
    NumOfComb = length(Comb);
    TempA = [];
    count = 0; %count number of variables
    for num = 1:NumOfComb
        i = Comb(num,1); j = Comb(num,2);
        EE = [E{i};E{j}]*U;
        for r1 = 1: height(EE)
            for r2 =1:height(EE)
                TempA = [TempA, vec(EE(r1,:)'*EE(r2,:))];
            end
        end
        count = count + height(EE)^2;
        K_In.s = [K_In.s height(EE)];
    end
    A_In = [A_In,TempA];
        
    b_In = vec(C);
    c_In=[-b;zeros(count,1)];%maximize
        
    %Outer Approximation
    %Solve in Sedumi dual form
    b_Out = -vec(C);%minimize
    c_Out = b; %linear constraints %will add more later on
    K_Out.f = m;
    K_Out.s=[];
    
    %Outer-Linear constraints
    At_Out = [];
    for i = 1:m% m is number of constraints
        At_Out = [At_Out; vec(A{i})'];
    end
    
    %Outer-symmetric constraints
    temp = ones(dim_mat);
    temp = tril(temp,-1);
    idx = find(temp); %index of symmetric element
    for i = 1:length(idx)
        temp_mat = zeros(dim_mat);
        temp_mat(idx(i)) = 1;
        temp_mat = temp_mat + (-temp_mat');
        At_Out = [At_Out; vec(temp_mat)'];
    end
    c_Out = [c_Out;zeros(length(idx),1)];%Symmetric Constraints
    
    K_Out.f = m+length(idx);%Linear Constraints and Symmetric Constraints
    
    %Outer-Cone constraints
    %In each iteration, the only one change happens here
    TempAt=[];
    count = 0; %count number of variables
    for num = 1:NumOfComb
        i = Comb(num,1); j = Comb(num,2);
        EE = [E{i};E{j}]*U;
        tempAt = [];
        for r1 = 1: width(EE)
            for r2 =1:width(EE)
                tempAt = [tempAt, -vec(EE(:,r1)*EE(:,r2)')];
            end
        end
        count = count + height(EE)^2;
        K_Out.s = [K_Out.s height(EE)];
        TempAt=[TempAt;tempAt];
    end    
    At_Out = [At_Out;TempAt];
    At_Out = At_Out';
    c_Out = [c_Out;zeros(count,1)];
    
    Iter = 0;
    OBJ_Out = []; %keep track of Obj Value
    OBJ_In = [];
    OBJ_U = [];
    while true
        if Iter > 0
            %change right hand side 
            TempA = [];
            count = 0; %count number of variables
            for num = 1:NumOfComb
                i = Comb(num,1); j = Comb(num,2);
                EE = [E{i};E{j}]*U;
                for r1 = 1: height(EE)
                    for r2 =1:height(EE)
                        TempA = [TempA, vec(EE(r1,:)'*EE(r2,:))];
                    end
                end
                count = count + height(EE)^2;
            end
            A_In(:,m+1:end) = [TempA];
        end
        
        [x_In,y_In,info_In]=sedumi(A_In,b_In,c_In,K_In);
        OBJ_In = [OBJ_In,-c_In(1:m)'*x_In(1:m)];
        
        
        [x_u,y_u,info_u,OBJ_u]=BFW2(A,b,C,P,U);
        OBJ_U = [OBJ_U,OBJ_u];
        


        %Outer Approximation
        if Iter>0
            TempAt=[];
            for num = 1:NumOfComb
                i = Comb(num,1); j = Comb(num,2);
                EE = [E{i};E{j}]*U;
                tempAt = [];
                for r1 = 1: width(EE)
                    for r2 =1:width(EE)
                        tempAt = [tempAt, -vec(EE(:,r1)*EE(:,r2)')];
                    end
                end
                TempAt=[TempAt;tempAt];
            end   
            At_Out(:,m+length(idx)+1:end) = TempAt';
%             At_Out = [At_Out;TempAt];
%             At_Out = At_Out';
        end
        
        %Solve in dual form
       [x_Out,y_Out,info_Out]=sedumi(At_Out,b_Out,c_Out,K_Out);
        
        X_Out = mat(y_Out);
        OBJ_Out = [OBJ_Out, trace(C*X_Out)];%keep track of OBJ value
        D = eig(X_Out);
        if Mode == 1
            if (min(D)>=0 || abs(min(D))<epsilon) 
                PSD = true; break;
            end
        elseif Mode == 2
            gap = abs(OBJ_U(end)-OBJ_Out(end))
            if (gap<epsilon || Iter>=20)  
                break;
            end
        end
        Iter = Iter +1;
        
        
        %find a new rotation
        S = vec(C);
        for i = 1:m
            S = S - x_In(i)*vec(A{i});
        end
        S = mat(S);
        D = eig(S);
        %since here is inner approximation, the matrix is guaranteed to be PSD
        if (sum(D<0)>0) %address numerical issue 
            S = S + 0.00001*eye(dim_mat);
        end
%         S = full(S);
%  %      Spectral decomposition******
%         [V D]=eig(S);
%         Norm = sqrt(sum(V.*V));
%         V = V./(ones(dim_mat,1)*Norm);
%         U = (V*sqrt(D))';
        
         U = chol(S);%cholesky decomposition *****
    end
end
