% construct model, control cage and deformed cage
UpperCircle = @ (x) sqrt(1-x.^2);
LowerCircle = @ (x) -sqrt(1-x.^2);
x = [-1:0.001:1];
y_pos = UpperCircle(x);
y_neg = LowerCircle(fliplr(x));
V = [x(1:end-1),fliplr(x);y_pos(1:end-1),y_neg]';
VC = [1.5,0;1.5,1.5;0,1.5;-1.5,1.5;-1.5,0;-1.5,-1.5;0,-1.5;1.5,-1.5];
FC = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,1];
VCD = VC;
VCD(3,:) = [0,0.8];
VCD(4,:) = [-1,0.6];
FCD = FC;
%% plot original
createfigure(VC)
hold on; plot(V(:,1),V(:,2),"LineWidth",3,"Color","b")
%%
[phi,psi,N,UN] = green_coordinates(VC,FC,V);
[~,~,N_FCD,~] = green_coordinates(VCD,FCD,V);
U = phi*VCD + psi*N_FCD;
createfigure(VCD)
hold on; plot(U(:,1),U(:,2),"LineWidth",3,"Color","g")


%%
W = MVC(V, VC);
VD = W * VCD;
plot(VD(:,1),VD(:,2),"LineWidth",3,"Color","m");




