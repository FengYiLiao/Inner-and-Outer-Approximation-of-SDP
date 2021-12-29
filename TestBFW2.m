%Test Block-Factor-width-2
%Inner approximation
[At,b,c,X,A,C] = Generate_SDP_Problems(5,10);

dim_mat = 5; %dimension of the symmetric matrix

P=[1 1 2 1]; %partition. Sum(p) = dim_mat

[x_fw2,y_fw2,info_fw2,OBJ_fw2,X_fw2]=BFW2(A,b,C,P);

disp(X_fw2);

