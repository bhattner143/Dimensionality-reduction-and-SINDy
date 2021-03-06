clc
clear all
close all
%% Read data from the file generated by VICON
filename='trial1_all_markers_Trajectories_100.xlsx';
M = xlsread(filename);
Morigin=M(5:end,3:5);
M=M(5:end,12:end);
[mM,nM]=size(M);
ViconFrameRate=100;
tdata=0:0.01:(mM-1)/ViconFrameRate';
numOfaxes=3;
%create a structure containing the unfiltered and filtered postion of 2nd column markers
f1='MarkerTimedata';
f2='MarkerMovement';
f3='MarkerMovSize';
f4='MarkerOriginMovData';
markers=struct(f1,tdata',f4,Morigin,f3,[mM nM],f2,M);
markerfields=fieldnames(markers);
%% Fill missing cells by interpolation
 for jj=1:nM
      jj
    markerObj=iddata((markers(1).(markerfields{4})(:,jj)),[],0.01);
     markerObj1=misdata(markerObj);
     markerObj1=detrend(markerObj1);
     markerobj2=fft(markerObj1);
     temp1(:,jj)= markerObj1.OutputData;
     temp2(:,jj)=markerobj2.OutputData;
 end
 markers(1).MarkerDetrendedMovData=temp1;
  markers(1).MarkerFFTMovData=temp2;
% save('eg2a_missingdata_interpolation_workspace2');
plot(1:100,abs(markers(1).MarkerFFTMovData(1:100,50)))