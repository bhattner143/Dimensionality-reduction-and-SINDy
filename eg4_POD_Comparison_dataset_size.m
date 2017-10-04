clc
clear all
close all
%% Load
load('eg3_preprocesing_workspace2');
%%
dataSelStart=500%436;
dataSelEnd=1000:500:2500;%636;
Y=cell(size(dataSelEnd,2),1);
V=cell(size(dataSelEnd,2),1);
for jjj=1:size(dataSelEnd,2)
    NumCoor=2;
    NumMark=size(markerfields,1);
    MarkersZ=markers.MarkerFilteredMovData(dataSelStart:dataSelEnd(jjj),:);
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
    %% SVD
    [U, Sigma, V{jjj}] = svd(Xaug, 'econ');
    sizeV(jjj)=(size(V{jjj}(1,:),2));
    Y{jjj}=U.'*Xaug;
     
     figure(1)
     hold on
     ax1=subplot(1,2,1);plot(1:size(V{1}(1,:),2),V{jjj}(1:size(V{1}(1,:),2),1),'Linewidth',2);
     
     hold off
     hold on
     ax2=subplot(1,2,2);plot(1:size(V{1}(2,:),2),V{jjj}(1:size(V{1}(1,:),2),2),'Linewidth',2);
     hold off
end
set(gca,'linewidth',1.5,...
     'xcolor',[0,0,0]);
setFigProp2([24.4,12],11);
I=legend(ax1 ,{strcat('${\bf v}_{1}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(1))),...
    strcat('${\bf v}_{1}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(2))),...
    strcat('${\bf v}_{1}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(3))),...
    strcat('${\bf v}_{1}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(4))),...
     },'FontSize',9);
     set(I,'interpreter','latex');
 xlabel(ax1,'$ t $','interpreter','latex');    
     I=legend(ax2 ,{strcat('${\bf v}_{2}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(1))),...
    strcat('${\bf v}_{2}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(2))),...
    strcat('${\bf v}_{2}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(3))),...
    strcat('${\bf v}_{2}(t)$ for ${\bf X}_{aug}\in {4020\times }$',num2str(sizeV(4))),...
     },'FontSize',9);
     set(I,'interpreter','latex');
      xlabel(ax2,'$ t $','interpreter','latex');
     %matlabToLatexEps('f7',300);
% figure
%     subplot(1,2,1),plot(1:size(Y{1:jjj}(1:2,:),2),Y{jjj}(1,:),'Linewidth',2);
%     subplot(1,2,2),plot(1:size(Y{1:jjj}(1:2,:),2),Y{jjj}(2,:),'Linewidth',2);
