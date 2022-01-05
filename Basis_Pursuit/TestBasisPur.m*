% This file domenstrates how to use Basis Pursuit of both Inner and Outer approximation.
% The below script demonstrates two examples.
% The algorithm solve in the form of 
% min tr(C*X)
% s.t. tr(A_i,X) = b_i, i = 1,2,...,m
%     X is PSD
%Main functions: DD_BasisPur_Outer, SDP, and DD_BasisPur_Inner. 
clc;clear;

%Start with a simple example. 2 Linear constraints. 3x3 design matrix.  
C = [2 1 0;1 2 1;0 1 2];
b = [1 ;0.5];
A = cell(2,1); %Linear Constraints are stored in a cell.
A{1} = eye(3);
A{2} = ones(3);

[x_Out,y_Out,info,OBJ_Out] = DD_BasisPur_Outer(A,b,C);
X_Out = mat(y_Out);%notice capital letter refer to matrix 
disp('X_Out'); disp(X_Out);

%Compare with solution solved directly by PSD cone
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
X_Sdp = mat(x_sdp);%notice capital letter refer to matrix 
disp('X_Sdp'); disp(X_Sdp);

%plot
plot(OBJ_SDP-OBJ_Out,'-o');
xlabel('Iteration');
ylabel('SDP - DD BasisPursuit Outer Approximation');
legend('DD');
grid on;


%% A general example
%This example takes a lot of time to compute
%Generate a problem of a 5x5 design matrix and 10 linear constraints 
[At,b,c,X,A,C] = Generate_SDP_Problems(5,10);
% We only need A, C, and b

[x_Out,y_Out,info,OBJ_Out] = DD_BasisPur_Outer(A,b,C);
X_Out = mat(y_Out);%notice capital letter refer to matrix 
disp('X_Out'); disp(X_Out);

%Compare with solution solved directly by PSD cone
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
X_Sdp = mat(x_sdp);%notice capital letter refer to matrix 
disp('X_Sdp'); disp(X_Sdp);

%plot
plot(OBJ_SDP-OBJ_Out,'-o');
xlabel('Iteration');
ylabel('SDP - DD BasisPursuit Outer Approximation');
legend('DD');
grid on;

%%
%Start with a simple example. 2 Linear constraints. 3x3 design matrix.  
C = [2 1 0;1 2 1;0 1 2];
b = [1 ;0.5];
A = cell(2,1); %Linear Constraints are stored in a cell.
A{1} = eye(3);
A{2} = ones(3);

[x_In,y_In,info_In,OBJ_In] = DD_BasisPur_Inner(A,b,C);

%Compare with solution solved directly by PSD cone
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
X_Sdp = mat(x_sdp);%notice capital letter refer to matrix 
disp('X_Sdp'); disp(X_Sdp);

%plot
plot(OBJ_In-OBJ_SDP,'-o');
xlabel('Iteration');
ylabel('SDP - DD BasisPursuit Inner Approximation');
legend('DD');
grid on;

%% A general example
%This example takes a lot of time to compute
%Generate a problem of a 5x5 design matrix and 10 linear constraints 
[At,b,c,X,A,C] = Generate_S
DP_Problems(5,10);
% We only need A, C, and b

[x_In,y_In,info_In,OBJ_In] = DD_BasisPur_Inner(A,b,C);

%Compare with solution solved directly by PSD cone
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
X_Sdp = mat(x_sdp);%notice capital letter refer to matrix 
disp('X_Sdp'); disp(X_Sdp);

%plot
plot(OBJ_In-OBJ_SDP,'-o');
xlabel('Iteration');
ylabel('SDP - DD BasisPursuit Inner Approximation');
legend('DD');
grid on;



