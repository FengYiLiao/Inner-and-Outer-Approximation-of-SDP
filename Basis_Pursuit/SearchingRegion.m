%Generate searching region
%  I + xA + yB, where A and B are 10x10 symmetric matrices
%  In each iteration, 
%  max x+y
%  s.t. I + xA + yB is PSD/FW_{a,2}^{n}/FW_{a,2}^{n}(U_{k})
%  U_{k} = chol(I+x^{*}A+y{*}B)

%Goal variables: PSD, BFW

%randomly symmetric matrics
rng(1);
n = 10;%dimension of decision matrix
Struc= ones(n);
A = sprandsym(Struc);
B = sprandsym(Struc);
x = linspace(-2*10^(-1),2*10^(-1),100); %adjust the interval to get finer graph
y = x ;
PSD = zeros(length(x));
BFW = cell(1,3);
P = [2,2,2,2,2];
%construct Operator (substract sub-principle matrix)
% combination: 
CardOfP = length(P);
Comb = nchoosek(1:CardOfP,2);%block factor-width-2
NumOfComb = length(Comb);
I = eye(n);
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

U = eye(n);%first U

for iter = 1:3 %run 3 times
    BFW{iter} =zeros(length(x));
    %Set up the obj, constranits, and cone for checking if the (x,y) is in the current cone
    [AA,cc,kk]=setup(E,Comb,A,B,U);
    i = 1 ; j = 1;
    for b = y
       for a = x
           %only check psd region in first iteration
           if (iter ==1)
               S = I + a*A+b*B;
               D = eig(S);
               if (min(D)>=0)
                   PSD(j,i) = 1;
               end
           end
           
           [x_bfw,y_bfw,info_bfw] = CheckBFW2(A,B,AA,cc,kk,a,b);
           if (info_bfw.pinf == 0)
               BFW{iter}(j,i) = 1;
           end
           i = i+1;
       end
       j = j+1;
       i = 1;
    end
    [x_bfw,y_bfw,info_bfw,U]=UpdateU(E,Comb,A,B,U); 
end

%plot searching region
FeasibleRegion(PSD,BFW);



function [AA,cc,K] = setup(E,Comb,A,B,U)
%Input:
%      E: a cell containing all the operators
%      Comb: combination
%      A,B: 10x10 symmetric matric
%      U: basis pursuit matrix
%
%The problem will be set up in the form:
%X = [a b c d .........]. Put all elements from each small PSD matric together. 
%AA = [vec((u1)(u1)') ...... vec((uk)(uk)')]]
%AA constains all the outer products.
%c: does not matter, we are not solving the optimization problem. The goal is to check if current (x,y) is in the current cone. 
%Example) cone: X is 4x4. Partition P = [1,2,1].
%               X = [e1,e2,e3][a b c [e1,e2,e3]' + [e1,e4][j k [e1,e4]' 
%                              d e f                       l m]  
%                              g h i]  
%                   +[e2,e3,e4][n o p  [e2,e3,e4]'
%                               q r s
%                               t u v]
%                 = a(e1)(e1)'+ b(e1)(e2)'+.......+v(e4)(e4)' 
    NumOfComb = length(Comb);
    n = width(A);
    count = 0;%count number of variables
    AA = [];
    K.s = [];
    for num = 1:NumOfComb
        i = Comb(num,1); j = Comb(num,2);
        EE = [E{i};E{j}]*U;
        for r1 = 1: height(EE)
            for r2 =1:height(EE)
                AA = [AA, vec(EE(r1,:)'*EE(r2,:))];
            end
        end
        count = count + height(EE)^2;
        K.s = [K.s height(EE)];
    end
    cc = ones(count,1);
end

function [X,Y,info] = CheckBFW2(A,B,AA,cc,kk,x,y)
    n = width(A);
    I = eye(n);
    S = I+x*A+y*B;
    bb = vec(S);
    %solve in primal form
    [X,Y,info] = sedumi(AA,bb,cc,kk);
end

function [x,y,info,NewU] = UpdateU(E,Comb,A,B,U)
%solve in sedumi dual form
%Decision variables:[x,y,a,b,c,......]
%At=[vec(A) vec(B) vec((u1)(u1)') ...... 
%    0,0,           I(NumOfVariables in sub-PSD matrices) ]'
%c= [vec(I(size of A))
%      0(NumOfVariables in sub-PSD matrices)    ]
%b=[1,1,0,0,0,.....]
%Given P=[a1,a2,a3,...]
%k.f=SizeOfA k.s=[a1+a2,a2+a3,....]: each element is the dimension of
%sub-PSD matrix.
    NumOfComb = length(Comb);
    n = width(A);
    At = [];
    k.s = [];
    count = 0; %number of variables
    for num = 1:NumOfComb
        i = Comb(num,1); j = Comb(num,2);
        EE = [E{i};E{j}]*U;
        for r1 = 1: height(EE)
            for r2 =1:height(EE)
                At = [At, vec(EE(r1,:)'*EE(r2,:))];
            end
        end
        count = count + height(EE)^2;
        k.s = [k.s , height(EE)];
    end
    At = [-vec(A),-vec(B),At];
    At = [At;zeros(count,2),-eye(count)];
    At = At';
    c = [vec(eye(n));zeros(count,1)];
    b = zeros(2+count,1);
    b(1) = 1; b(2) = 1;
    k.f = n^2; 
    [x,y,info] = sedumi(At,b,c,k);
    I = eye(n);
    S = I+y(1)*A+y(2)*B;
    NewU = chol(S);
end

function FeasibleRegion(PSD,BFW)
%PSD is a matrix
%BFW is a cell containing a sequence of matrices
h = height(PSD);
w = width(PSD);
C = ones(h,w,3);
C1 = ones(h,w);C2= ones(h,w);C3 = ones(h,w);
C1(PSD>0) = 0;
C2(PSD>0) = 1;
C3(PSD>0) = 0;

C1(BFW{1}>0) = 0;
C2(BFW{1}>0) = 0;
C3(BFW{1}>0) = 1;

C1(BFW{2}>0) = 1;
C2(BFW{2}>0) = 0;
C3(BFW{2}>0) = 1;

C1(BFW{3}>0) = 0;
C2(BFW{3}>0) = 0;
C3(BFW{3}>0) = 0;

C(:,:,1) = C1;C(:,:,2) = C2;C(:,:,3) = C3;
image(C);
end
