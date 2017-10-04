clc
clear all
close all
%% LOAD MARKER DATA AND REDUCED NUMBER OF STATES OBTAINED FROM POD
load('eg9_SVD_Marker_datav2')
timedata=markers.MarkerTimedata(1:1451);
lambda =0.1; 
Sizelambda=size(lambda,2);
nState=2; 
FinPt=550;   
MSErr=zeros(Sizelambda,2);
Xistore=struct();
%%
 for iii=1:Sizelambda 
    %% INITIALISE TIME AND STATES
    xt=Vr(1:end-1,:);
    timedata1=timedata(1:FinPt);
    dt=timedata1(2)-timedata1(1);
    [m,n]=size(xt);% m=time stamp %n number of states
    %% INITIALIZE REGRESSION PARAMETERS
       % lambda is our sparsification knob.
    polyorder =2;  % search space up to fifth order polynomials
    usesine = 0;    % no trig functions
    %% COMPUTE DERIVATIVE
    time_diff = 'FD'
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
    Theta = poolData(xt(1:FinPt,:),n,polyorder,usesine);
    m3 = size(Theta,2);
    %% compute Sparse regression: sequential least squares
    Xi = sparsifyDynamics(Theta,dx(1:FinPt,:),lambda(iii),n);
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
    %% Plotting for the following parameter value
 figure(1)
    kkk=1.5;
    timedata_xB=tB;
    subplot(2,1,1),plot(timedata(1:end-1),xt(:,1),'Linewidth',kkk);
    hold on;plot(timedata_xB,xB(:,1),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex');
    text(0,0.055,strcat('$\lambda$=',num2str(lambda(iii)),', Training dataset size=',num2str(FinPt),', Poly order=',num2str(polyorder)),'interpreter','latex','FontSize',12)
    setFigProp2([15,15],12)
    set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
    I=legend({strcat('$x_1$'),...
        strcat('$\tilde {x}_1 $')},...
        'FontSize',12,'Location','northeast');
    set(I,'interpreter','latex')
    subplot(2,1,2),plot(timedata(1:end-1),xt(:,2),'Linewidth',kkk);
    setFigProp2([15,15],12)
    hold on;plot(timedata_xB,xB(:,2),'Linewidth',kkk);xlabel('time ($s$)','interpreter','latex')
%     text(12.5,10.8,strcat('$\lambda$=',num2str(lambda(iii))),'interpreter','latex','FontSize',9)
    I=legend({strcat('$x_2$'),...
        strcat('$\tilde {x}_2 $')},...
        'FontSize',12,'Location','northeast');
    set(I,'interpreter','latex')
    set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
 hold off
% if lambda(iii)==0.1||lambda(iii)==0.15||lambda(iii)==0.2
%   filename = sprintf('%s_%d','PDE_find_withot_RIP',iii)
   matlabToLatexEps('Prediction_polyorder2c',300);
% end
 end
%% NRMSE vs lambda
figure
semilogy(lambda,MSErr,'Linewidth',2)
xh=xlabel('$\lambda$','interpreter','latex');
% set(gca, 'XTick', [0.01:0.05:0.3]);
yh=ylabel('$\log_{10} NRMSE$','interpreter','latex');
I=legend({'$NRMSE_{(x_1^{actual},x_1^{predicted})}$','$NRMSE_{(x_2^{actual},x_2^{predicted})}$'},'FontSize',12);
set(I,'interpreter','latex')
% set(gca, 'XTick', 0.1:0.3:3);
set(gca,'linewidth',1,...
     'xcolor',[0,0,0]);
%set(gca,'XTickLabel',sprintf('%0.1f\n',lambda));
%matlabToLatexEps('NRMSE_lambda0',300);
% set(gca,'YTickLabel',sprintf('%1.4f\n',err));