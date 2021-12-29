%%SOCP using Sedumi
%%dual standard form
clc;

% %example 1
% dim_mat = 3; %dimension of the symmetric matrix
% b = -[2 2 2 0 2 2]';
% dim_Var = length(b);
% At = -[1 0 1 0 0 1;1 2 1 2 2 1]';%constraint
% c = [-1,-0.5]';
% K.f = 2;
% K.l = 0;
% K.q =[];
% OBJ =[];

% %example 2
% dim_mat = 3; %dimension of the symmetric matrix
% b = -[1 4 9 6 0 7]';
% dim_Var = length(b);
% At = -[1 0 3 2 14 5;0 4 6 16 0 4]';%constraint
% c = [-11,-19]';
% K.f = 2;
% K.l = 0;
% K.q =[];
% OBJ =[];



%example 3
dim_mat = 2; %dimension of the symmetric matrix
b = -[2 2 0]';
dim_Var = length(b);
At = -[1 0 1]';%constraint
c = [-1]';
K.f = 1;
K.l = 0;
K.q =[];
OBJ =[];

Comb = nchoosek(1:dim_mat,2);
NumOfComb = length(Comb);
if dim_mat <3
    NumOfComb = 1;
end
for j = 1:NumOfComb
    V = zeros(dim_mat,2);
    V(Comb(j,1),1)=1;
    V(Comb(j,2),2)=1;
    [C1 C2 C3] = ConvMat2SOCPCons(V,dim_Var,dim_mat);
    At = [At(:,1:K.f+K.l),-C1',-C3',At(:,K.f+K.l+1:end),-(C1+C3)',-2*C2',-(C1-C3)'];
    c = [c;zeros(5,1)];
    K.l= K.l+ 2;
    K.q =[K.q 3];
end

idx = Subscribpts(dim_mat);
PSD = false;
Iter = 1;
err = 1.0e-3;
while ~PSD
    % Perform the optimization.
    [x,y,info]=sedumi(At,b,c,K);
    OBJ = [OBJ, b'*y]; 
    X = zeros(dim_mat);
    X(idx) = y;% Show the optimal x solution.
    X = X + tril(X,-1)'
    [V,D] = eig(X);
    [d, I]=sort(diag(D));
    if d(1)>=0  %get PSD 
        PSD = true;
        break;
    end
    v = V(:,I(1:2));
    v(:,1)=v(:,1)/norm(v(:,1));
    v(:,2)=v(:,2)/norm(v(:,2));
    Iter = Iter+1;
    [C1 C2 C3] = ConvMat2SOCPCons(v,dim_Var,dim_mat);
    At = [At(:,1:K.f+K.l),-C1',-C3',At(:,K.f+K.l+1:end),-(C1+C3)',-(C1-C3)',-2*C2'];
    c = [c;zeros(5,1)];
    K.l= K.l+ 2;
    K.q =[K.q 3];%cautious
end
disp("Iter")
disp(Iter)
Y=zeros(dim_mat);
Y(idx) = y;
X = Y + tril(Y,-1)';
plot(OBJ);