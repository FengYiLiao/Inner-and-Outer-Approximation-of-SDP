%Basis pursuit LP
%Outer approximation
%Basis Pursuit Outer Approximation relies on the inner approximation of
%generating a sequence of rotation matrices, we denote it as U in the
%script. 

function  [x_Out,y_Out,info_Out,OBJ]=DD_BasisPur_Outer(A,b,C)
%Input:
%      A: a cell storing a sequence of symmetric matrices.
%      b: a column vector.
%      C: the coefficient matrix of objective function

    epsilon = 10^(-3);
    dim_mat = width(C); %dimension of C;
    m = length(A);%number of linear constraints
    z = 1; %index of Extreme Ray
    v = zeros(dim_mat,dim_mat^2);%Original extreme Raies (n^2)
    
    %generate extreme ray (1)
    Comb1 = nchoosek(1:dim_mat,1);
    NumOfComb1 = length(Comb1);
    for j = 1:NumOfComb1
        u = zeros(1,dim_mat)';
        u(j) = 1;
        v(:,z) = u; z = z +1;
    end

    %generate extreme ray (2)
    Comb2 = nchoosek(1:dim_mat,2);
    NumOfComb2 = length(Comb2);
    if dim_mat <3 % can be deleted
        NumOfComb2 = 1;
    end
    
    for j = 1:NumOfComb2
        for i = 1:2
            u = zeros(1,dim_mat)';
            if (i == 1)
                u(Comb2(j,1))=1;
                u(Comb2(j,2))=-1;
            else
                u(Comb2(j,1))=1;
                u(Comb2(j,2))=1;
            end
            v(:,z) = u; z = z +1;
        end
    end
    
    U = eye(dim_mat); %the first U

    
    PSD = false;
    
    %start with inner approximation
    %Solve in Sedumi dual form
    b_In =  [b; zeros(dim_mat^2,1)]; %decision variables: y and alphas
    c_In = [vec(C);zeros(dim_mat^2,1)];%Linear constraints and nonnegative constraints
    
    %Cone constraints
    At_In = [];
    for i = 1:m% m is number of constraints
        At_In = [At_In, vec(A{i})];
    end
    
    %Cone constraints- Right hand side
    %In each iteration, the only one change happens here
    TempAt = [];
    for i = 1:width(v)
        TempAt = [TempAt, vec((U'*v(:,i))*((U'*v(:,i))'))];
    end
    At_In = [At_In, TempAt];
    K_In.f = dim_mat^2;
    
    %nonnegative constraints
    TempAt= [zeros(dim_mat^2,m), -eye(dim_mat^2)]; %y and alphas
    At_In = [At_In; TempAt];
    At_In = At_In';
    K_In.l = dim_mat^2;
    
    %Outer Approximation
    %Solve in Sedumi dual form
    b_Out = -vec(C);
    c_Out = [b ;zeros(dim_mat^2,1)]; %linear constraints and nonnegative constraint
    
    %Outer-Linear constraints
    At_Out = [];
    for j = 1:m% m is number of constraints
        At_Out = [At_Out, vec(A{j})];
    end
    %Outer-symmetric constraints
    temp = ones(dim_mat);
    temp = tril(temp,-1);
    idx = find(temp); %index of symmetric element
    for i = 1:length(idx)
        temp_mat = zeros(dim_mat);
        temp_mat(idx(i)) = 1;
        temp_mat = temp_mat + (-temp_mat');
        At_Out = [At_Out, vec(temp_mat)];
    end
    c_Out = [c_Out;zeros(length(idx),1)];%Symmetric Constraints
    
    K_Out.f = m+length(idx);%Linear Constraints and Symmetric Constraints
   
    %Outer-Cone constraints
    %In each iteration, the only one change happens here
    temp = [];
    for i = 1:width(v)
        temp = [temp, -vec((U'*v(:,i))*(U'*v(:,i))')];
    end
    At_Out = [At_Out,temp];
    K_Out.l = dim_mat^2;
        
    Iter = 0;
    OBJ = [];
    while ~PSD
        if Iter > 0
            %change right hand side 
            TempAt = [];
            for i = 1:width(v)
                TempAt = [TempAt, vec((U'*v(:,i))*((U'*v(:,i))'))];
            end
            TempAt = [TempAt;-eye(dim_mat^2)];
            %this is the only change
            TPS_At_In = At_In';
            TPS_At_In(:,m+1:m+dim_mat^2) = TempAt;
            At_In = TPS_At_In';
        end
        %notice that we didn't impose symmetric constraints for Inner Approximation
        %but the way we construct construct the problem had already imposed it.
        [x_In,y_In,info_In]=sedumi(At_In,b_In,c_In,K_In);

        %find a new rotation
        S = vec(C);
        for i = 1:m
            S = S - y_In(i)*vec(A{i});
        end
        S = mat(S);
        D = eig(S);
        %since here is inner approximation, the matrix is guaranteed to be PSD
        if (sum(D<0)>0) %address numerical issue 
            S = S + 0.00001*eye(dim_mat);
        end

        U = chol(S);%cholesky decomposition *****

        if Iter > 0
            %Change Cone constraints
            temp = [];
            for i = 1:width(v)
                temp = [temp, -vec((U'*v(:,i))*(U'*v(:,i))')];
            end
            At_Out(:,m+length(idx)+1:m+length(idx)+dim_mat^2) = temp;        
        end

        [x_Out,y_Out,info_Out]=sedumi(At_Out,b_Out,c_Out,K_Out);
        X_Out = mat(y_Out);
        OBJ = [OBJ, trace(C*X_Out)];%keep track of OBJ value
        D = eig(X_Out);
        if (min(D)>=0 || abs(min(D))<epsilon) 
            PSD = true; break;
        end
        Iter = Iter +1;
    end
end


