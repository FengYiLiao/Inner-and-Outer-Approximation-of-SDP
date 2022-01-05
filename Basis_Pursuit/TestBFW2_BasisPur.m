%Test BFW_BasisPur_Inner
%Start with a simple example. 2 Linear constraints. 3x3 design matrix.  
% C = [2 1 0;1 2 1;0 1 2];
% b = [1 ;0.5];
% A = cell(2,1); %Linear Constraints are stored in a cell.
% A{1} = eye(3);
% A{2} = ones(3);
% P = [1,2];
format long
k = 10;%dimensions
m = 10;%number of linear constraints
[At,b,c,X,A,C] = Generate_SDP_Problems(k,m);
P1 = [2,2,2,2,2];
[x_In_BFW2_1,y_In_BFW2_1,info,OBJ_In_BFW2_1] = BFW2_BasisPur_Inner(A,b,C,P1);
P2 = [4,2,2,2];
[x_In_BFW2_2,y_In_BFW2_2,info,OBJ_In_BFW2_2] = BFW2_BasisPur_Inner(A,b,C,P2);
P3 = [4,4,2];
[x_In_BFW2_3,y_In_BFW2_3,info,OBJ_In_BFW2_3] = BFW2_BasisPur_Inner(A,b,C,P3);

[x_In_SDD,y_In_SDD,info,OBJ_In_SDD] = BFW2_BasisPur_Inner(A,b,C,ones(1,k));
[x_In_DD,y_In_DD,info_In_DD,OBJ_In_DD] = DD_BasisPur_Inner(A,b,C);
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);
% X_Out = mat(y_Out);%notice capital letter refer to matrix 
% disp('X_Out'); disp(X_Out);

%plot
plot(OBJ_In_BFW2_1-OBJ_SDP,'-o');
hold on;
plot(OBJ_In_BFW2_2-OBJ_SDP,'-o');
plot(OBJ_In_BFW2_3-OBJ_SDP,'-o');
plot(OBJ_In_SDD-OBJ_SDP,'-o');
plot(OBJ_In_DD-OBJ_SDP,'-o');
hold off;



xlabel('Iteration');
ylabel('SDP-BWF2/SDD/DD');
legend('BFW2_{P1}','BFW2_{P2}','BFW2_{P3}','SDD','DD');
t =['n=' ' ' num2str(k) ', m=' ' ' num2str(m)];
title(t);
grid on;
