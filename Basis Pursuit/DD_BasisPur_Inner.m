%Basis pursuit LP
%Inner approximation

function  [x,y,info,OBJ]=DD_BasisPur_Inner(A,b,C)
%Input:
%A: a cell storing a sequence of symmetric matrices.
%b: a column vector.
%C: the coefficient matrix of objective function
%Sovle in Sedumi dual form
%y=[x_11,x_12,...,x_nn,a_1,a_2,...,a_n^2]
    dim_mat = width(C);%dimension of design matirx
    m = length(A);%number of linear constraints

    z = 1; %index of Extreme Ray
    v = zeros(dim_mat,dim_mat^2);%Original extreme Ray n^2
    %generate extreme ray (1)
    Comb1 = nchoosek(1:dim_mat,1);
    NumOfComb1 = length(Comb1);
    for j = 1:NumOfComb1
        u = zeros(1,dim_mat)';
        u(j) = 1;
        v(:,z) = u;z = z +1; 
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

    b_new = [-vec(C); zeros(dim_mat^2,1) ];%X and alphas
    c_new = zeros(m+2*dim_mat^2,1);%m linear constraints
    %dim_mat cone constraints. dim_mat nonnegative constraints
    c_new(1:m) = b;%linear constraints

    %linear constraints
    At_new = [];
    for j = 1:m% m is number of constraints
        temp = [vec(A{j});zeros(dim_mat^2,1)];%since we have 2*dim_mat design variables
        At_new = [At_new, temp];
    end

    %Cone constraints
    TempAt = -eye(dim_mat^2);
    for i = 1:width(v)
        %TempAt = [TempAt, vec(ERay{i})] ;
        TempAt = [TempAt, vec((U'*v(:,i))*(U'*v(:,i))')] ;
    end
    TempAt= TempAt';
    At_new = [At_new, TempAt];
    K_new.f = m+dim_mat^2;

    %non-negative constraints
    TempAt= [zeros(dim_mat^2), -eye(dim_mat^2)]';
    At_new = [At_new, TempAt];
    K_new.l = dim_mat^2;

    [x,y,info]=sedumi(At_new,b_new,c_new,K_new); %The first iteration
    X = y(1:dim_mat^2);%y=[x_11,x_12,...,x_nn,a_1,a_2,...,a_n^2] we only need the first n^2
    X = mat(X);
    CNVG = false;%Convergence or Not
    %Convergence condition happens when Norm(abs(PreObj-NowObj)) < err
    Iter = 1; err = 1.0e-3;
    OBJ = []; OBJ = [OBJ, trace(C*X)];%keep track of Obj Value
    PreX = X ; 
    epsilon = 0.00001;
    while ~CNVG
        D = eig(PreX);
        %since here is inner approximation, the matrix is guaranteed to be PSD
        if (sum(D<0)>0)%address numerical issue
            PreX = PreX + epsilon*eye(dim_mat);
        end
        
        U = chol(PreX);% naming follows the paper
        
        %Cone constraints
        %In each iteration, the only change happens here
        TempAt = -eye(dim_mat^2);
        for i = 1:width(v)
            %TempAt = [TempAt, vec(ERay{i})] ;
            TempAt = [TempAt, vec((U'*v(:,i))*(U'*v(:,i))')] ;
        end
        TempAt= TempAt';
        At_new(:,m+1:m+dim_mat^2) = TempAt; %
        
        [x,y,info]=sedumi(At_new,b_new,c_new,K_new);
        X = y(1:dim_mat^2);
        X = mat(X);
        OBJ = [OBJ trace(C*X)];
        %if (trace(eye(dim_mat)*abs((PreX-X)))) < err
        if  (OBJ(Iter)-OBJ(Iter+1))^2<err
            CNVG = true;
            break;
        end
        PreX = X;
        Iter = Iter+1;
        
    end 
end


