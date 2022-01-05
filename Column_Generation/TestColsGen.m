% This file demonstrates how to use LPColsGen and SOCPColsGen
% The below script demonstrates two examples.
% Both algorithms are doing outer approximation of PSD Cone
% The algorithm solve in the form of 
% min C'*X
% s.t. A*X = b
% Note: 
%     X is cast as a vector constaing only a half part of the elements in the
%     matrix 
%     A is cast as a matrix containing all the linear constraints
%     This format follows the startd LP
%Main functions: LPColsGen, SDP, and SOCPColsGen. 


%% Start with a simple example
clear,clc;

%for LP and SOCP
c = [2 2 2 0 2 2]';
b = [1,0.5]';
At = [1 0 1 0 0 1;1 2 1 2 2 1];%constraint
K.f = 2; %the cone should be specified in advance

%for SDP
A = cell(2,1); %Linear Constraints are stored in a cell.
A{1} = eye(3);
A{2} = ones(3);
C = [2 1 0;1 2 1;0 1 2];

[x,y,info,OBJ_LP]=LPColsGen(At,b,c,K);%%Compute Optimization
[x_SOCP,y_SOCP,info_SOCP,OBJ_SOCP]=SOCPColsGen(At,b,c,K);%%Compute Optimization

%show the decision matrix %%LP
[M N]=size(At);
dim_mat = (-1+sqrt(1+8*N))/2; %dimension of the symmetric matrix
idx = Subscribpts(dim_mat);%get the indices of the desgin elements in the matrix
Y=zeros(dim_mat);
Y(idx) = y;
X = Y + tril(Y,-1)'

%Compare with SDP
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
plot(OBJ_SDP-OBJ_LP,'-o');
hold on ;
plot(OBJ_SDP-OBJ_SOCP,'-o');
xlabel('Iteration');
ylabel('SDP - LP/SOCP');
legend('LP','SOCP');
grid on;
hold off;


%% A general example
m = 10;%number of linear constraints
k = 5;%dimension of design matrix
[At,b,c,X,A,C]=Generate_SDP_Problems(k,m);

K.f = m; %the cone should be specified in advance

[x,y,info,OBJ_LP]=LPColsGen(At,b,c,K);%%Compute Optimization
[x_SOCP,y_SOCP,info_SOCP,OBJ_SOCP]=SOCPColsGen(At,b,c,K);%%Compute Optimization

%show the decision matrix %%LP
[M N]=size(At);
dim_mat = (-1+sqrt(1+8*N))/2; %dimension of the symmetric matrix
idx = Subscribpts(dim_mat);%get the indices of the desgin elements in the matrix
Y=zeros(dim_mat);
Y(idx) = y;
X = Y + tril(Y,-1)'

%Compare with SDP
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
plot(OBJ_SDP-OBJ_LP,'-o');
hold on ;
plot(OBJ_SDP-OBJ_SOCP,'-o');
xlabel('Iteration');
ylabel('SDP - LP/SOCP');
legend('LP','SOCP');
grid on;
hold off;

