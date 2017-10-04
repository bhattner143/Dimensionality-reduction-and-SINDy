
clc
clear all
close all
numcyc=2;
T=2.25;
A=100;
tf=numcyc*T;
t = 0:0.01:tf;
y =A* sin(2*pi*t / T);

m=1;
n=1;

switch m
    case 1 
        %% Generate data at 200ms
        
        for ii=1:round(0.2/0.01):round(tf/0.01)
            if y(ii)>=0
                y1(n)=y(ii);
                t1(n)=t(ii);
            else
                y1(n)=0;
                t1(n)=t(ii);
            end
        n=n+1;
        end
%% Generate data at 500ms
    case 2
        
        for ii=1:round(0.5/0.01):round(tf/0.01)
         if y(ii)>=0
                y1(n)=y(ii);
                t1(n)=t(ii);
            else
                y1(n)=0;
                t1(n)=t(ii);
            end
         n=n+1;
        end
end