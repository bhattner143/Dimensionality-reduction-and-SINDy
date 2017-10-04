clc
clear all
close all
%% LOAD MARKER DATA AND REDUCED NUMBER OF STATES OBTAINED FROM POD
load('eg9_SVD_Marker_datav2')
timedata=markers.MarkerTimedata(1:1451);
lambda =0.01:0.01:0.5; %Polyorder2 0.01:0.01:0.5 Polyorder3 0.01:0.01:3
Sizelambda=size(lambda,2);
nState=2; 
FinPt=400:50:900;   %Polyorder2 400:50:900 Polyorder3 200:50:1000
SizeDataset=size(FinPt,2);
MSErr=zeros(SizeDataset,Sizelambda,2);
Xistore=struct();
for jjj=1:SizeDataset
%%
jjj
 for iii=1:Sizelambda 
     iii
    %% INITIALISE TIME AND STATES
    xt=Vr(1:end-1,:);
    timedata1=timedata(1:FinPt(jjj));
    dt=timedata1(2)-timedata1(1);
    [m,n]=size(xt);% m=time stamp %n number of states
    %% INITIALIZE REGRESSION PARAMETERS
       % lambda is our sparsification knob.
    polyorder =2;  % search space up to fifth order polynomials
    usesine = 0;    % no trig functions
    %% COMPUTE DERIVATIVE
    time_diff = 'FD';
     if strcmp(time_diff,'poly' )==1
            m2 = m-2*width_t;
            offset_t = width_t;
     else  
            n2 = n;
            offset_t = 0;
     end
    m2=m;
    dx=zeros(m2+1,n2,'double');
    kk=1;
        for i=kk:nState
            dx(:,i)=FiniteDiff2(Vr(:,i),dt,1);
        end
    dx=dx(1:end-1,:);
    %% pool Data  (i.e., build library of nonlinear time series)
    Theta = poolData(xt(1:FinPt(jjj),:),n,polyorder,usesine);
    m3 = size(Theta,2);
    %% compute Sparse regression: sequential least squares
    Xi = sparsifyDynamics(Theta,dx(1:FinPt(jjj),:),lambda(iii),n);
    NumSparseCoeff(jjj,iii)=nnz(~Xi);
    fieldname1=strcat('Eoptimised',num2str(FinPt(jjj)));
    fieldname2=strcat('NumSparseCoeff',num2str(FinPt(jjj)));
    Xistore(iii).(fieldname1)=Xi;
    Xistore(iii).(fieldname2)=nnz(~Xi);
    %% Prediction
    xt0=xt(1,:)';
    options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n)); 
    x=xt;
    tspan=timedata(1:end-1);
    temp=zeros(size(tspan,1),3);
    [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,xt0,options);  % approximate
    MaxErr=max(MSErr);
    if size(xB(:,1),1)==size(xt(:,1),1)
        Err = abs(xt-xB).^2;
        MSErr(jjj,iii,:) = sqrt(sum(Err))/size(xt,1);
    else
        MSErr(jjj,iii,:) =[MaxErr(1,1) MaxErr(1,1)];
    end
    %% Plotting for the following parameter value
%  figure
%     kkk=1.5;
%     title(strcat('Plots for polyorder=',num2str(polyorder),'and use of sine term=',num2str(usesine),'value of lambda',num2str(lambda(iii))));
%     timedata_xB=tB;
%     subplot(2,1,1),plot(timedata(1:end-1),xt(:,1),'Linewidth',kkk);
%     hold on;plot(timedata_xB,xB(:,1),'Linewidth',kkk);xlabel(strcat('$\lambda$=',num2str(lambda(iii)),' ,FinPt=',num2str(FinPt(jjj))),'interpreter','latex')%xlabel('time ($s$)','interpreter','latex');
%     text(1.5,1.5,strcat('$\lambda$=',num2str(lambda(iii))),'interpreter','latex','FontSize',9)
%     setFigProp2([15,15],11)
%     set(gca,'linewidth',1,...
%      'xcolor',[0,0,0]);
%     I=legend({'$x_1^{actual}$','$x_1^{predicted}$ '},'FontSize',9,'Location','northeast');
%     set(I,'interpreter','latex')
%     subplot(2,1,2),plot(timedata(1:end-1),xt(:,2),'Linewidth',kkk);
%     setFigProp2([15,15],11)
%     hold on;plot(timedata_xB,xB(:,2),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex')
%     text(1.5,1.5,strcat('$\lambda$=',num2str(lambda(iii))),'interpreter','latex','FontSize',9)
%     I=legend({'$x_2^{actual}$','$x_2^{predicted}$'},'FontSize',9,'Location','northeast');
%      set(I,'interpreter','latex')
%     set(gca,'linewidth',1,...
%      'xcolor',[0,0,0]);
%     subplot(3,1,3),plot(timedata,xt1(:,3),'Linewidth',kkk);
%     setFigProp2([15,15],11)
%     hold on;plot(timedata_xB,xB(:,3),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex')
%      text(12.5,10.8,strcat('$\lambda$=',num2str(lambda(iii))),'interpreter','latex','FontSize',9) 
%     I=legend({'$x_3^{actual}$','$x_3^{predicted}$'},'FontSize',9,'Location','northeast');
%      set(I,'interpreter','latex')
% set(gca,'linewidth',1,'xcolor',[0,0,0]);
% if lambda(iii)==0.1||lambda(iii)==0.15||lambda(iii)==0.2
%   filename = sprintf('%s_%d','PDE_find_withot_RIP',iii)
%   matlabToLatexEps(filename,300);
% end
 end
end
%% NRMSE vs lambda vs datset size
CM=jet(SizeDataset);
for jjj=1:SizeDataset
    figure(1)
   
  semilogy(lambda,(MSErr(jjj,:,1)),'color',rand(1,3),'Linewidth',1.5); 
    hold on
   legendInfo{jjj} = ['Datset size = ' num2str(FinPt(jjj))];
end
I=legend(legendInfo,'FontSize',12,'Location','eastoutside')
%  set(I,'interpreter','latex')
xh=xlabel('$\lambda$','interpreter','latex');
yh=ylabel('$NRMSE{(\bf {x}_1,\bf {\tilde x}_1)}$','interpreter','latex');
set(gca,'linewidth',1,'xcolor',[0,0,0]);
setFigProp2([15,15],12)
hold off
  filename = 'NRMSE_Polyorder2'
  matlabToLatexEps(filename,300);
[xx,yy]=meshgrid(lambda,FinPt);
figure(2)
subplot(1,2,1),surf(xx,yy,MSErr(:,:,1));%colormap winter
xh=xlabel('$\lambda$','interpreter','latex');
yh=ylabel('Size of the dataset','interpreter','latex');
zh=zlabel('$NRMSE{(x_1^{actual},x_1^{predicted})}$','interpreter','latex');
% I=legend({'$NRMSE_{(x_1^{actual},x_1^{predicted})}$','$NRMSE_{(x_2^{actual},x_2^{predicted})}$'},'FontSize',12);
% set(I,'interpreter','latex')
% set(gca, 'XTick', 0.1:0.3:3);
set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
 hold on
 surf(xx,yy,NumSparseCoeff(:,:,1));
 hold off
 subplot(1,2,2),surf(xx,yy,MSErr(:,:,2));%colormap winter
xh=xlabel('$\lambda$','interpreter','latex');
yh=ylabel('Size of the dataset','interpreter','latex');
zh=zlabel('$NRMSE{(\bf {x}_1,x_1^{predicted})}$','interpreter','latex');
% I=legend({'$NRMSE_{(x_1^{actual},x_1^{predicted})}$','$NRMSE_{(x_2^{actual},x_2^{predicted})}$'},'FontSize',12);
% set(I,'interpreter','latex')
% set(gca, 'XTick', 0.1:0.3:3);
set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
 % %matlabToLatexEps('NRMSE_lambda0',300);