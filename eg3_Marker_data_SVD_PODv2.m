
clc
clear all
close all
%% Load
load('eg3_preprocesing_workspace2');
%%
dataSelStart=500%436;
dataSelEnd=2000%636;
NumCoor=2;
NumMark=size(markerfields,1);
MarkersZ=markers.MarkerFilteredMovData(dataSelStart:dataSelEnd,:);
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
%% Economic SVD
[U, Sigma, V] = svd(Xaug, 'econ');
%% r rank truncation
r=2; % rank truncation
for iii=1
    iii
Ur = U(:, 1:r(iii)); 
Sigma_r = Sigma(1:r(iii), 1:r(iii));
Vr = V(:, 1:r(iii));
Xaug_r=Ur*Sigma_r*Vr';
ProjXaugOnPODmode1=U(:,1).'*Xaug;% also equal to 1st SV times 1st column of V 
%where each column of V represents the time dynamics corresponding to POD mode 1
ProjXaugOnPODmode2=U(:,2).'*Xaug;% also equal to 2nd SV times 2st column of V
ProjXaugOnPODmode3=U(:,3).'*Xaug;
ProjXaugOnPODmode4=U(:,4).'*Xaug;
Y=U.'*Xaug;
%Error
temp=sum(sum((Xaug-Xaug_r).^2));
Err(iii)=temp./(size(Xaug,1)*size(Xaug,2));
end
% save('eg9_SVD_Marker_datav2');
%% Effect od dataset size on reduced time dynamics

%% Plotting
close all
figure (1)
surfl(abs(X(:,:)));shading interp; colormap winter; colorbar
view(33,60)
xlabel('Time snapshots, $m$','interpreter','latex')
ylabel('States, $x$','interpreter','latex')
zlabel('${\bf x}(t)/\bf X$','interpreter','latex')
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0]);
setFigProp2([12.2,6],11);

%   matlabToLatexEps('f1',300);
figure (2)
surfl(abs(Xaug(:,:))); shading interp; colormap winter; colorbar
view(33,60)
xlabel('Time snapshots, $m_{aug}$','interpreter','latex')
ylabel('States, $n_{aug}$','interpreter','latex')
zlabel('$|{\bf X}_{aug}|$','interpreter','latex')
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0],...
          'fontsize',11);
setFigProp2([24.4,12],11);
% matlabToLatexEps('f2',600);
figure(3)
plot(diag(Sigma)/sum(diag(Sigma)),'ko','Linewidth',2);
ylabel('Normalised SV, diag($\Sigma$)/sum(diag($\Sigma$))','interpreter','latex')
xlabel('POD Mode ')
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0],...
          'fontsize',9);
setFigProp2([12.2,6],9); 
%  matlabToLatexEps('f3',600);
figure (4)
plot(abs(U(1:120,1:2)),'Linewidth',2);% plotting the first two modes
I=legend('POD mode 1 (${\bf u}_{1})$','POD mode 2 (${\bf u}_{2})$');
set(I,'interpreter','latex');
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0],...
          'fontsize',11);
setFigProp2([12.2,6],11);      
% matlabToLatexEps('f4',600);
figure(5)
surfl(abs(Xaug_r(:,:)));shading interp; colormap winter; colorbar
view(33,60)
xlabel('Snapshots, $m_{aug}$','interpreter','latex')
ylabel('States, $n_{aug}$','interpreter','latex')
zlabel('${\bf X}_{aug}^{r}$','interpreter','latex')
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0],...
          'fontsize',11);
setFigProp2([24.4,12],11);
% matlabToLatexEps('f5',600);
figure(6);
surfl(Y); shading interp; colormap winter; colorbar %how the modes are projected in the new basis coordinate system
xlabel('$m_{aug}/t$','interpreter','latex')
ylabel('POD modes,${\bf u}_{i}$','interpreter','latex')
zlabel('$\bf V$','interpreter','latex')
view(-78,12)
setFigProp2([24.4,12],11);         
% matlabToLatexEps('f6',300);
% how the first four modes are projected in the new basis coordinate system
figure(7);
nn=jj;
plot(1:nn,ProjXaugOnPODmode1,1:nn,ProjXaugOnPODmode2,1:nn,ProjXaugOnPODmode3,1:nn,ProjXaugOnPODmode4,'Linewidth',1);
I=legend({'Time dynamics (${\bf v}_{1}(t)$) corresponding to POD mode 1 (${\bf u}_{1}$)',...
    'Time dynamics  (${\bf v}_{2}(t)$) corresponding to POD mode 2 (${\bf u}_{2}$)',...
    'Time dynamics  (${\bf v}_{3}(t)$) corresponding to POD mode 3 (${\bf u}_{3}$)',...
    'Time dynamics  (${\bf v}_{4}(t)$) corresponding to POD mode 4 (${\bf u}_{4}$)'},...
    'FontSize',5,'Location','eastoutside');
xlabel('$t$','interpreter','latex')
set(I,'interpreter','latex');
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0],...
          'fontsize',11);
setFigProp2([12.2,6],11);         
% matlabToLatexEps('f7a',300);
figure(8)
subplot(3,1,1),surf(X(:,1:1451))
;shading interp; colormap winter; colorbar
view(33,60)
subplot(3,1,2),surf(Xaug(1:60,:))
;shading interp; colormap winter; colorbar
view(33,60)
subplot(3,1,3),surf(Xaug_r(1:60,1:1451))
;shading interp; colormap winter; colorbar
view(33,60)