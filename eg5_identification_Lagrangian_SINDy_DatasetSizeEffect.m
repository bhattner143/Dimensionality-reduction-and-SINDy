clc
clear all
close all
%% LOAD SAVED PROCESSED MARKER DATA
load('eg9_SVD_Marker_datav2')
timedata=markers.MarkerTimedata(1:1451);
lambda =0.1; 
nState=2;
FinPt=100:50:750;   
SizeDataset=size(FinPt,2);
MSErr=zeros(SizeDataset,2);
Xistore=struct();
%%
 for iii=1:SizeDataset
    %% More processing
    xt=Vr(1:end-1,:);
    timedata1=timedata(1:FinPt(iii));
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
    Theta = poolData(xt(1:FinPt(iii),:),n,polyorder,usesine);
    m3 = size(Theta,2);
    %% compute Sparse regression: sequential least squares
    Xi = sparsifyDynamics(Theta,dx(1:FinPt(iii),:),lambda,n);
    Xistore(iii).Eoptimised=Xi;
    %% Prediction
    xt0=xt(1,:)';
    options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n)); 
    x=xt;
    tspan=timedata(1:end-1);
    temp=zeros(size(tspan,1),3);
    [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,xt0,options);  % approximate
    if size(xB(:,1),1)==size(xt(:,1),1)
        Err = abs(xt-xB).^2;
        MSErr(iii,:) = sqrt(sum(Err))/size(xt,1);
    else
        MSErr(iii,:) =[1 1];
    end
    %% Plotting
    % kk=1;
    % kk1=1.5
    % for i=1:nMar
    %     figure
    %     for i=kk:kk-1+size(CoorSel,2)
    %     subplot(1,1,i-kk+1),plot(timedata,M12Processed(:,i),'r',timedata,dx(:,i),'Linewidth',kk1);
    %     xlabel('time in seconds');
    %     legend(strcat('x_',num2str(i),'(t)'),strcat('xdot_',num2str(i),'(t)'));
    %     end
    %     kk=kk+size(CoorSel,2);
    % end
    %  plot(M12Processed(:,1),dx(:,1));xlabel('x_1');ylabel('xdot_1');
    % hold on
    % plot(M12Processed(:,2),dx(:,2));
    % hold off
    %  plot(tspan,x(:,2));hold on
    %   plot(tspan,xB(:,2))
    %% Plotting for the following parameter value
 figure
    kkk=1.5;
    title(strcat('Plots for polyorder=',num2str(polyorder),'and use of sine term=',num2str(usesine),'value of lambda',num2str(lambda)));
    timedata_xB=tB;
    subplot(2,1,1),plot(timedata(1:end-1),xt(:,1),'Linewidth',kkk);
    hold on;plot(timedata_xB,xB(:,1),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex');
    text(12.5,2.3,strcat('$\lambda$=',num2str(FinPt(iii))),'interpreter','latex','FontSize',9)
    setFigProp2([15,15],11)
    set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
    I=legend({'$x_1^{actual}$','$x_1^{predicted}$ '},'FontSize',9,'Location','northeast');
    set(I,'interpreter','latex')
    subplot(2,1,2),plot(timedata(1:end-1),xt(:,2),'Linewidth',kkk);
    setFigProp2([15,15],11)
    hold on;plot(timedata_xB,xB(:,2),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex')
    text(12.5,10.8,strcat('$\lambda$=',num2str(FinPt(iii))),'interpreter','latex','FontSize',9)
    I=legend({'$x_2^{actual}$','$x_2^{predicted}$'},'FontSize',9,'Location','northeast');
     set(I,'interpreter','latex')
    set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
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
figure
semilogy(1:SizeDataset,MSErr)
% %% NRMSE vs lambda
% figure
% for iii=1:3
%     hold on
%     plot(lambda,err(:,iii),'Linewidth',1.5);
%        setFigProp2([15,8],11)
% end
% xh=xlabel('$\lambda$','interpreter','latex');
% yh=ylabel('$NRMSE$ without RIP','interpreter','latex');
% I=legend({'$NRMSE_{(x_1^{actual},x_1^{predicted})}$','$NRMSE_{(x_2^{actual},x_2^{predicted})}$','$NRMSE_{(x_3^{actual},x_3^{predicted})}$'},'FontSize',12);
%   set(I,'interpreter','latex')
% set(gca, 'XTick', [0.01:0.05:0.3]);
% set(gca,'linewidth',1,...
%      'xcolor',[0,0,0]);
% %matlabToLatexEps('NRMSE_lambda0',300);
% % set(gca,'XTickLabel',sprintf('%1.1f\n',lambda));
% % set(gca,'YTickLabel',sprintf('%1.4f\n',err));
