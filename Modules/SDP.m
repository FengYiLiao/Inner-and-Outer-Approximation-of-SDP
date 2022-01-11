function [x,y,info,OBJ] = SDP(A,b,C)
%Solve in Sedumi "Primal" form
%Input:
%      A: a cell storing a sequence of symmetric matrices.
%      b: a column vector.
%      C: the coefficient matrix of objective function
    c_new = vec(C);
    b_new = b;
    dim_mat = width(C);
    A_new = zeros(length(A),length(c_new));
    for j = 1: length(A)
       A_new(j,:) = vec(A{j})'; 
    end
    K.s = [dim_mat];
    [x,y,info] = sedumi(A_new,b_new,c_new,K);
    X = mat(x);
    OBJ = trace(C'*X);
end