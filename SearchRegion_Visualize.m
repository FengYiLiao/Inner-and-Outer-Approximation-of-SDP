%Improve visualization of searching region
%Please load the matrices first (search region)
load('SearchRegion_2_New_BFW2.mat'); %Load the corresponding matrices
fill_clr = [0.803921580314636   0.878431379795074   0.968627452850342
            0.992156863212585   0.917647063732147   0.796078443527222
            0.839215695858002   0.909803926944733   0.850980401039124
            0.937254905700684   0.866666674613953   0.866666674613953];
line_clr = [0.000000000000000   0.447058826684952   0.741176486015320;
            0.850980401039124   0.325490206480026   0.098039217293262;
            0.000000000000000   0.498039215803146   0.000000000000000;
            0.494117647409439   0.184313729405403   0.556862771511078];
[y,x] = find(PSD);
k = convhull(x,y);
%xt = x*0.008+(-2*10^(-1));
%yt = y*0.008+(-2*10^(-1));
xt = x*0.012+(-3*10^(-1));
yt = y*0.012+(-3*10^(-1));
%plot(xt(k),yt(k),'k');
patch(xt(k),yt(k),'-', 'LineWidth', 1, 'FaceColor',fill_clr(1,:) , 'EdgeColor',line_clr(1,:));
hold on;
for n = 1:2
   [y,x] = find(BFW{n}); %Change 
   %[y,x] = find(DD{n});
   k = convhull(x,y);
   %plot(x(k),y(k));
   %xt = x*0.008+(-2*10^(-1));
   %yt = y*0.008+(-2*10^(-1));
   xt = x*0.012+(-3*10^(-1));
   yt = y*0.012+(-3*10^(-1));
   patch(xt(k),yt(k),'-', 'LineWidth', 1, 'FaceColor',fill_clr(n+1,:) , 'EdgeColor',line_clr(n+1,:));
   hold on 
end
hold off
xlabel('x')
ylabel('y')
alpha(.5)
xlim([-0.3,0.35])
ylim([-0.3,0.4])
%grid on
%xlim([-0.2,0.3])
%ylim([-0.2,0.3])
%legend('PSD','$DD$','$DD(U_1)$','Interpreter','latex','Location','north','EdgeColor','none','Orientation','horizontal')
%legend('PSD','$SDD$','$SDD(U_1)$','Interpreter','latex','Location','north','EdgeColor','none','Orientation','horizontal')
legend('PSD','$FW^{10}_{\alpha,2}$','$FW^{10}_{\alpha,2}(U_1)$','Interpreter','latex','Location','north','EdgeColor','none','Orientation','horizontal')