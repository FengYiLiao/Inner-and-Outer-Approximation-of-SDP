function [C1 C2 C3]=ConvMat2SOCPCons(V,n,k)
    %Input: n*2 vector V, n: dimension of x, k: dimension of x in matrix form
    %V has N*M dimension %In SOPCColumnGeneration m=2
    %Output: coefficient of the constraints
    %We compute the matrix   [C1     C2
    %                         C2     C3]
    %where C1, C2,and C3 are the coefficients of the variables  
    
    [N M] = size(V);
    G =[];%constraints matrix 
    for j=1:M
       m = [];
       for i=1:M
           test = V(:,j)*V(:,i)';
           m =[m ConvMat2LinerCons(test)'];
       end
       G=[G;m];
    end
    C1 = G(1,1:n);
    C2 = G(2,1:n);
    C3 = G(2,n+1:2*n);
end

