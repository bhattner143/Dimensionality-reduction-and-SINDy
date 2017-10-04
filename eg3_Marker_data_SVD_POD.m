clc
clear all
close all
%% Load
load('eg3_preprocesing_workspace');
%%
xdata=[1 2 3 4 5];
dataSelStart=500%436;
dataSelEnd=2000%636;
NumCoor=2;
NumMark=size(markerfields,1);
MarkersZ=zeros(dataSelEnd-dataSelStart+1,NumMark);
index1=1;
for j=1:NumMark
    MarkersZ(:,index1:index1+1)=markers(4).(markerfields{j})(dataSelStart:dataSelEnd,2:3);
    index1=index1+2;
end
X=MarkersZ.';
[nX, mX]=size(X);
%% hankel matrix
hX=2*size(X,2)/size(X,1);
hX=ceil(hX);
Xaug=zeros(hX*nX,mX-hX+1);
temp1=1;
for ii=1:hX
   
    temp2=nX*ii;
    for jj=1:mX-hX+1
        Xaug(temp1:temp2,jj)  =X(:,ii+jj-1);  
    end
    temp1=temp2+1;
end
%% SVD and rank-2 truncation
r =161; % rank truncation
[U, Sigma, V] = svd(Xaug, 'econ');
Ur = U(:, 1:r); 
Sigma_r = Sigma(1:r, 1:r);
Vr = V(:, 1:r);
Xaug_r=Ur*Sigma_r*Vr';
Y=pinv(U(:,:))*Xaug;
ProjXaugOnPODmode1=U(:,1).'*Xaug;% also equal to 1st SV times 1st column of V 
%where each column of V represents the time dynamics corresponding to POD mode 1
ProjXaugOnPODmode2=U(:,2).'*Xaug;% also equal to 2nd SV times 2st column of V
ProjXaugOnPODmode3=U(:,3).'*Xaug;
ProjXaugOnPODmode4=U(:,4).'*Xaug;
Y=U.'*Xaug;
Y2=Xaug(:,1)*U(:,1).';
%% Plotting
close all
figure (1)
surfl(abs(X(:,:)));
xlabel('Time snapshots, $m$','interpreter','latex')
ylabel('States, $n$','interpreter','latex')
zlabel('Amplitude','interpreter','latex')
figure (2)
surfl(abs(Xaug(:,:))); %shading interp
xlabel('Time snapshots, $m$','interpreter','latex')
ylabel('States, $n_{aug}$','interpreter','latex')
zlabel('Amplitude','interpreter','latex')
figure(3)
plot(diag(Sigma)/sum(diag(Sigma)),'ko','Linewidth',2);
ylabel('Normalised SV, diag($\Sigma$)/sum(diag($\Sigma$))','interpreter','latex')
xlabel('POD Mode')
figure (4)
plot(abs(U(1:30,1:3)),'Linewidth',2);% plotting the first three modes
I=legend('POD mode 1','POD mode 2','POD mode 3');
figure(5)
surfl(abs(Xaug_r(:,:)));% shading interp
xlabel('Snapshots, $m$','interpreter','latex')
ylabel('States, $n_aug$','interpreter','latex')
zlabel('Amplitude','interpreter','latex')
figure(6);
surfl(Y); %shading interp% how the modes are projected in the new basis coordinate system
xlabel('Time dynamics of the corrresponding POD mode, $m$','interpreter','latex')
ylabel('POD modes,${\bf u}_{i}$','interpreter','latex')
zlabel('Amplitude','interpreter','latex')
figure;
% how the first four modes are projected in the new basis coordinate system
nn=1201;
plot(1:nn,ProjXaugOnPODmode1,1:nn,ProjXaugOnPODmode2,1:nn,ProjXaugOnPODmode3,1:nn,ProjXaugOnPODmode4,'Linewidth',2);
I=legend('Time dynamics corresponding to POD mode 1',...
    'Time dynamics corresponding to POD mode 2',...
    'Time dynamics corresponding to POD mode 3',...
    'Time dynamics corresponding to POD mode 4');