%Test BFW2_BasisPur_Outer

k = 10;%dimensions
m = 10;%number of linear constraints
[At,b,c,X,A,C] = Generate_SDP_Problems(k,m);
P1 = ones(1,k);
Mode = 2;
[x_Out_BFW2,y_Out_BFW2,info,OBJ_Out_BFW2,OBJ_In_BFW2] = BFW2_BasisPur_Outer(A,b,C,P1,Mode);

X = reshape((y_Out_BFW2),[k,k])
plot(OBJ_In_BFW2,'-or');
hold on ;
plot(OBJ_Out_BFW2,'-ob');
plot(OBJ_In_BFW2-OBJ_Out_BFW2,'-o');
legend('In','Out','Gap'); hold off;
%%

%Compare with solution solved directly by PSD cone
[x_sdp,y_sdp,info_sdp,OBJ_SDP]=SDP(A,b,C);

%plot
plot(OBJ_SDP-OBJ_Out_BFW2,'-o');
xlabel('Iteration');
ylabel('SDP - DD BasisPursuit Outer Approximation');
legend('DD');
grid on;