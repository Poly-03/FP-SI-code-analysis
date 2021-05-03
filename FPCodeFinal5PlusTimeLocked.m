%%
filename='933 TEST BOTTOM2020-06-10T12_54_38.csv'; %input filename to be analyzed
keysfilename='933 TEST KEYSB2020-06-10T12_54_38.csv'; %input keyfilename to be analyzed
%%
%read in data

d=readtable(filename); %reads in csv
k=readtable(keysfilename);
d(:,4)=[]; %clip superfluous 4th column
k(:,3)=[]; %clip superfluous 3th column
d=table2array(d); %convert to an array

filename=filename(1:20); %clip filename extension   STOP CODE HERE AND CORRECT K table if corrections needed 
%%
    %find data column/time column in data table
    if mean(d(:,1))<1000000
        datad=d(:,1);
    end

    if mean(d(:,2))<1000000
        datad=d(:,2);
    end

    if mean(d(:,3))<1000000
        datad=d(:,3);
    end

    if mean(d(:,1))>1000000 && mean(d(:,1))<100000000
        time=d(:,1);
    end

    if mean(d(:,2))>1000000 && mean(d(:,2))<100000000
        time=d(:,2);
    end

    if mean(d(:,3))>1000000 && mean(d(:,3))<100000000
        time=d(:,3);
    end

%%  
%housekeeping
    
    if rem(length(datad),2) == 1 %check to see if data is of an even or odd length, correct to even if odd, for deinterleaving
        datad(1,:)=[];
        time(1,:)=[];
    end
    
    
    clear m;
        
    FP.timestamp=time; %store uncorrected timestamps
    FP.data(:,1)=time; %put data and time into a single array
    FP.data(:,2)=datad; 
    
%      FP.data(1:1632,:)=[]; %clip problematic data ADJUST THIS!!!!!!
    %%
    figure
    plot(FP.data(:,1),FP.data(:,2)) % plot unprocessed data
        mkdir (filename)
            set(gcf, 'Position', [400, 300, 1000, 600])
     figA=gcf;
     figAstr='/Unprocessed.jpg';
     figAbstr='/Unprocessed.fig';
     filefigAstr=strcat(filename,figAstr);
     filefigAbstr=strcat(filename,figAbstr);
    saveas(figA,filefigAstr,'jpeg'); %save figure
    saveas(figA,filefigAbstr,'fig');
    
    %%
  
    keytimes=k(:,2);
    keytimes=table2array(keytimes);
    keytimes(:,1) = (keytimes(:,1) - time(1,1))/(60*1000); %convert key times to relative time in minutes
    
%%    
%normalize time and deinterleave
FP.data(:,1) = (FP.data(:,1) - time(1,1))/(60*1000); %convert to relative time in minutes

FP.deinterleave(:,1) = downsample(FP.data(1:end,1),2); %deinterleave -- in this case, column 2 = calcium, column 3 - isosbestic
FP.deinterleave(:,3) = downsample(FP.data(1:end,2),2); %if this code won't run you may have an indivisible length of data, just open the data FP.data and clip a row off the start.
FP.deinterleave(:,2) = downsample(FP.data(2:end,2),2);

% if mean(FP.deinterleave(:,2)) < mean(FP.deinterleave(:,3)) %switch order if 410/470 switched
%     FP.deinterleave(:,1) = downsample(FP.data(1:end,1),2); %deinterleave -- in this case, column 2 = calcium, column 3 - isosbestic
%     FP.deinterleave(:,3) = downsample(FP.data(1:end,2),2); %if this code won't run you may have an indivisible length of data, just open the data FP.data and clip a row off the start.
%     FP.deinterleave(:,2) = downsample(FP.data(2:end,2),2); %NOTE: This code doesn't always work - if you see the 470 is on the bottom in the original signal, omit this code
% end


  %%
  %plot raw data  
%   
 raw_reference = FP.deinterleave(2:end,3)';
 raw_signal = FP.deinterleave(2:end,2)';
% 

figure
subplot(2,1,1)
plot(FP.deinterleave(:,1),FP.deinterleave(:,2),'m')
 str = ['470-(',filename];
    title(str) 
    xlabel('Time (minutes)')
    ylabel('Mean Pixel Value')
    xlim([0 5.2])
subplot(2,1,2)
plot(FP.deinterleave(:,1),FP.deinterleave(:,3),'b')
 str = ['410-(',filename];
     title(str) 
    xlabel('Time (minutes)')
    ylabel('Mean Pixel Value')
    xlim([0 5.2])

    set(gcf, 'Position', [400, 300, 1000, 600])
     fig1=gcf;
     fig1str='/fig1.jpg';
     fig1bstr='/fig1.fig';
     filefig1str=strcat(filename,fig1str);
     filefig1bstr=strcat(filename,fig1bstr);
    saveas(fig1,filefig1str,'jpeg'); %save figure
    saveas(fig1,filefig1bstr,'fig');
    
    %%
    %z-score and convert to dFF
    %usage: zdFF = get_zdFF(reference, signal, smooth_win, remove, lambda, itermax, order, wep, p)
    
    % Calculates z-score dF/F signal based on fiber photometry calcium-idependent 
% and calcium dependent signals.
%
% This program is a translation in MATLAB of the Python source code of
% get_zdFF.py
%
% Input 
%     reference: calcium-independent signal (usually 405-420 nm excitation)
%     signal: calcium-dependent signal (usually 465-490 nm excitation 
%              for green fluorescent proteins, or ~560 nm for red)
%     smooth_win: window for moving average smooth 
%     remove: the beginning of the traces with a steep slope one would like to remove
%   Inputs for airPLS:
%     lambda: lambda is an adjustable parameter, it can be adjusted by user. 
%             The larger lambda is, the smoother baseline will be 
%     itermax: maximum iteration times
%     order: an integer indicating the order of the difference of penalties
%     wep: weight exception proportion at both the start and end
%     p: asymmetry parameter for the start and end
%
%  Output
%     zdFF - z-score dF/F, vector
%    
%  Examples:
%     zdFF = get_zdFF(reference, signal);
%     zdFF = get_zdFF(reference, signal, 10, 200, 5e9, 50, 2, 0.5, 0.5);
%
%  Reference:
%     (1) Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry 
%         to Record Neural Activity in Freely Moving Animal. J. Vis. Exp. 
%         (152), e60278, doi:10.3791/60278 (2019)
%         https://www.jove.com/video/60278/multi-fiber-photometry-to-record-neural-activity-freely-moving
%
%  March 2020 Ekaterina Martianova ekaterina.martianova.1@ulaval.ca



%     zdFF = get_zdFF(raw_reference,raw_signal,20,30); %using a 20-frame smoothing window here (~1s), ~1.5s clip-out at beginning (for 40hz interleaved recordings)
   
timeCorrected = FP.deinterleave(:,1);
timeCorrected(1:30)=[]; %must match the clip out length, also on line 227, remove must match
%     
%     figure
% plot(timeCorrected,zdFF','k')
%  str = ['zscore-(',filename];
%     title(str) 
%     xlabel('Time (minutes)')
%     ylabel('z-score')
    


%%
%Step by Step
%smooth data

smooth_win = 20; %smooth over a 20-frame window (equals ~1s)
smooth_reference = movmean(raw_reference,smooth_win);
smooth_signal = movmean(raw_signal,smooth_win);

figure
subplot(2,1,1)
plot(smooth_signal,'m')
 str = ['Smooth 470-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
subplot(2,1,2)
plot(smooth_reference,'b')
 str = ['Smooth 410-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
%%
%Remove slope with AirPLS function

lambda = 5e11;
order = 2;
wep = 0.1;
p = 0.5;
itermax = 50;

[reference,base_r]= airPLS(smooth_reference,lambda,order,wep,p,itermax);
[signal,base_s]= airPLS(smooth_signal,lambda,order,wep,p,itermax);

figure
subplot(2,1,1)
plot(smooth_signal,'m')
xlim([0 6100]);
hold on
plot(base_s,'k')
 str = ['AirPLS 470-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
hold off
subplot(2,1,2)
plot(smooth_reference,'b')
xlim([0 6100]);
hold on
plot(base_r,'k')
 str = ['AirPLS 410-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
hold off

  set(gcf, 'Position', [400, 300, 1000, 600])
     figA=gcf;
     figAstr='/AirPLS.jpg';
     figAbstr='/AirPLS.fig';
     filefigAstr=strcat(filename,figAstr);
     filefigAbstr=strcat(filename,figAbstr);
    saveas(figA,filefigAstr,'jpeg'); %save figure
    saveas(figA,filefigAbstr,'fig');
    
    clear figA
    clear figAstr
    clear figAbstr
    clear filefigAstr
    clear fielfigAbstr
%%
%Remove start of recording
remove = 30; %here was using 30 frames = first 1.5 seconds
reference = reference(remove:end);
signal = signal(remove:end);

figure
subplot(2,1,1)
plot(signal,'m')
 str = ['Clip Start 470-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
subplot(2,1,2)
plot(reference,'b')
 str = ['Clip Start 410-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
%%
%z-score data
z_reference = (reference - median(reference)) / std(reference);
z_signal = (signal - median(signal)) / std(signal);

figure
subplot(2,1,1)
plot(z_signal,'m')
 str = ['Z-Score Smooth 470-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
subplot(2,1,2)
plot(z_reference,'b')
 str = ['Z-Score Smooth 410-(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
    
      set(gcf, 'Position', [400, 300, 1000, 600])
     fig2=gcf;
     fig2str='/fig2.jpg';
     fig2bstr='/fig2.fig';
     filefig2str=strcat(filename,fig2str);
     filefig2bstr=strcat(filename,fig2bstr);
    saveas(fig2,filefig2str,'jpeg'); %save figure
    saveas(fig2,filefig2bstr,'fig');
%%
%Fit reference signal to calcium signal using non negative robust linear regression
fitdata = fit(z_reference',z_signal',fittype('poly1'),'Robust','on');

figure
hold on
plot(z_reference,z_signal,'k.')
plot(fitdata,'b')
 str = ['Fitted 470 + 410 -(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
hold off
%%

%Align Reference (410) to Signal (470)

z_reference = fitdata(z_reference)';

figure
plot(z_signal,'m')
 hold on
plot(z_reference,'b')
xlim([0 6100]);
 str = ['Aligned -(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
hold off

  set(gcf, 'Position', [400, 300, 1000, 600])
     fig3=gcf;
     fig3str='/fig3.jpg';
     fig3bstr='/fig3.fig';
     filefig3str=strcat(filename,fig3str);
     filefig3bstr=strcat(filename,fig3bstr);
    saveas(fig3,filefig3str,'jpeg'); %save figure
    saveas(fig3,filefig3bstr,'fig');
%%
%Calculate Z-Score dFF

FP.zdFF = z_signal - z_reference;
figure
plot(FP.zdFF,'k')
xlim([0 6100]);
 str = ['Z-dFF -(',filename];
    title(str) 
    xlabel('Time (frames)')
    ylabel('Mean Pixel Value')
   %% 
   %plot trends
   
    figure
    FP.zdFF2=FP.zdFF';
     FP.smoothed = smooth(FP.zdFF2(:,1),100);
       plot(timeCorrected,FP.smoothed, 'm', 'linewidth', 2);
       hold on
       plot(timeCorrected,FP.zdFF2, 'k');
 str = ['zscore-(',filename];
    title(str) 
    xlabel('Time (minutes)')
    ylabel('z-score')
    
    
    set(gcf, 'Position', [100, 400, 1800, 250])
      fig4=gcf;
     fig4str='/fig4.jpg';
     fig4bstr='/fig4.fig';
     filefig4str=strcat(filename,fig4str);
     filefig4bstr=strcat(filename,fig4bstr);
    saveas(fig4,filefig4str,'jpeg'); %save figure
    saveas(fig4,filefig4bstr,'fig');

    %%
    %House Keeping
    clear fig1
    clear fig1str
    clear fig1bstr
    clear filefig1str
    clear filefig1bstr
    clear fig2
    clear fig2str
    clear fig2bstr
    clear filefig2str
    clear filefig2bstr
    clear fig3
    clear fig3str
    clear fig3bstr
    clear filefig3str
    clear fileifg3bstr
    clear fig4
    clear fig4str
    clear fig4bstr
    clear filefig4str
    clear filefig4bstr
    
    %%
    %Sort Keystrokes
    
    k2=k(:,1);
    k2=table2array(k2);
    index = strfind(k2, 'S');
    keyindex(:,1) = cellfun('isempty',index);
    index = strfind(k2, 'D');
    keyindex(:,2) = cellfun('isempty',index);
    index = strfind(k2, 'F');
    keyindex(:,3) = cellfun('isempty',index);
    keyindex=not(keyindex);
    keyindex(end+1,:)=0;
    
    
    keydiffs=diff(keyindex);
    dindex=(find(keydiffs(:,2)==1));
    findex=(find(keydiffs(:,3)==1));
    dindex=dindex+1;
    findex=findex+1;
    dindex(:,2)=(find(keydiffs(:,2)==-1));
    findex(:,2)=(find(keydiffs(:,3)==-1));
    
         dtimelist=0;
    
      for i=1:length(dindex(:,1))
          number=num2str(i);
           name=strcat('dtimelistt_',number);
          dtimelistt(i,:)=dindex(i,1):dindex(i,2);
          str=[name,'(:,1)=dtimelistt(i,:);'];
          eval(str)
          clear dtimelistt
          clear str
          str=['dtimelist(end+1:(end+length(',name,')),',number,')=(',name,');'];
              eval(str)
              str=['clear ',name,';'];
              eval(str)
                                  
      end
      
      dtimelist=nonzeros(dtimelist);
              
               
    ftimelist=0;
    
      for i=1:length(findex(:,1))
          number=num2str(i);
           name=strcat('ftimelistt_',number);
          ftimelistt(i,:)=findex(i,1):findex(i,2);
          str=[name,'(:,1)=ftimelistt(i,:);'];
          eval(str)
          clear ftimelistt
          clear str
          str=['ftimelist(end+1:(end+length(',name,')),',number,')=(',name,');'];
              eval(str)
                 str=['clear ',name,';'];
              eval(str)
                    
      end
      
      ftimelist=nonzeros(ftimelist);
    
    dtimes=keytimes(dtimelist,1);
    ftimes=keytimes(ftimelist,1);
    
    FP.dtimes(:,1)=dtimes(1:2:end,1);
    FP.dtimes(:,2)=dtimes(2:2:end,1);
    
    FP.ftimes(:,1)=ftimes(1:2:end,1);
    FP.ftimes(:,2)=ftimes(2:2:end,1);
    
     FP.dzdFF=[0 0];
     FP.dzdFF2=[0 0];

%%
%Find segments that are <5s apart (0.0833 min) 
for i=2:length(FP.dtimes)
    FP.D_TimeCheck(i-1,1)=FP.dtimes(i,1);
    FP.D_TimeCheck(i-1,2)=FP.dtimes(i-1,2);
end
FP.D_TimeCheck(:,3)=FP.D_TimeCheck(:,1)-FP.D_TimeCheck(:,2);
FP.D_Segs=find((FP.D_TimeCheck(:,3))>0.0833);
FP.D_SegsMerge(1,1)=FP.dtimes(1,1);

for i=1:length(FP.D_Segs)
    FP.D_SegsMerge(i,2)=FP.D_TimeCheck(FP.D_Segs(i),2);
    FP.D_SegsMerge(i+1,1)=FP.D_TimeCheck(FP.D_Segs(i),1);
end
FP.D_SegsMerge(end,2)=FP.dtimes(end,2);

for i=2:length(FP.ftimes)
    FP.F_TimeCheck(i-1,1)=FP.ftimes(i,1);
    FP.F_TimeCheck(i-1,2)=FP.ftimes(i-1,2);
end

FP.F_TimeCheck(:,3)=FP.F_TimeCheck(:,1)-FP.F_TimeCheck(:,2);
FP.F_Segs=find((FP.F_TimeCheck(:,3))>0.0833);
FP.F_SegsMerge(1,1)=FP.ftimes(1,1);

for i=1:length(FP.F_Segs)
    FP.F_SegsMerge(i,2)=FP.F_TimeCheck(FP.F_Segs(i),2);
    FP.F_SegsMerge(i+1,1)=FP.F_TimeCheck(FP.F_Segs(i),1);
end
FP.F_SegsMerge(end,2)=FP.ftimes(end,2);

        %%
    % Add 2 seconds (0.033 min) to beginning of each interaction epoch,  &to
    % end of each epoch - use this for non timelocked analysis
     %FP.D_SegsMerge(:,1)=FP.D_SegsMerge(:,1)-0.033;
  %   FP.D_SegsMerge(:,2)=FP.D_SegsMerge(:,2)+0.033;
     
     %FP.F_SegsMerge(:,1)=FP.F_SegsMerge(:,1)-0.033;
   %   FP.F_SegsMerge(:,2)=FP.F_SegsMerge(:,2)+0.033;

  %%   
 for i=1:length(FP.D_SegsMerge)
     
      FP.dSignal=find(FP.D_SegsMerge(i,2) > timeCorrected(:,1) & timeCorrected(:,1) > FP.D_SegsMerge(i,1));
      number=num2str(i);
      name=strcat('FP.D.dzdFF_',number);
      str=[name,'(1:length(FP.dSignal),1)=(timeCorrected(FP.dSignal(1):FP.dSignal(end)));'];
      eval(str);
      clear str;
      name=strcat('FP.dzdFF2_',number);
      str=[name,'(1:length(FP.dSignal),1)=(FP.zdFF2(FP.dSignal(1):FP.dSignal(end)));'];
      eval(str);
      clear str;
      name=strcat('FP.D.dzdFF_',number);
      str=['FP.dzdFF(end+1:(end+length(',name,')),',number,')=(',name,');'];
      eval(str);
      clear str;
      
      name=strcat('FP.dzdFF2_',number);     
      str=['FP.dzdFF2(end+1:(end+length(',name,')),',number,')=(',name,');'];
      eval (str)
      
      clear str
      clear FP.dSignal;
    
 end
     FP.dzdFF=nonzeros(FP.dzdFF);
     FP.dzdFF2=nonzeros(FP.dzdFF2);
FP.D=[];
  
FP.D_Signal(:,2)=FP.dzdFF2;
FP.D_Signal(:,1)=FP.dzdFF;
      figure
      
plot(FP.D_Signal(:,1),FP.D_Signal(:,2),'c','linewidth', 2)
xlim([0 5.2]);
hold on
%%
     FP.fzdFF=[0 0];
     FP.fzdFF2=[0 0];
     
 for i=1:length(FP.F_SegsMerge)
     
      FP.fSignal=find(FP.F_SegsMerge(i,2) > timeCorrected(:,1) & timeCorrected(:,1) > FP.F_SegsMerge(i,1));
      number=num2str(i);
      name=strcat('FP.F.dzdFF_',number);
      str=[name,'(1:length(FP.fSignal),1)=(timeCorrected(FP.fSignal(1):FP.fSignal(end)));'];
      eval(str);
      clear str;
      name=strcat('FP.fzdFF2_',number);
      str=[name,'(1:length(FP.fSignal),1)=(FP.zdFF2(FP.fSignal(1):FP.fSignal(end)));'];
      eval(str);
      clear str;
      name=strcat('FP.F.dzdFF_',number);
      str=['FP.fzdFF(end+1:(end+length(',name,')),',number,')=(',name,');'];
      eval(str);
      clear str;
      
      name=strcat('FP.fzdFF2_',number);     
      str=['FP.fzdFF2(end+1:(end+length(',name,')),',number,')=(',name,');'];
      eval (str)
      
      clear str
      clear FP.fSignal;
    
 end
     FP.fzdFF=nonzeros(FP.fzdFF);
     FP.fzdFF2=nonzeros(FP.fzdFF2);
FP.F=[];
  
FP.F_Signal(:,2)=FP.fzdFF2;
FP.F_Signal(:,1)=FP.fzdFF;



plot(FP.F_Signal(:,1),FP.F_Signal(:,2), 'g', 'linewidth', 2)

figure
plot(timeCorrected(:,1),FP.zdFF2(:,1),'k');
xlim([0 5.2]);
hold on
plot(FP.F_Signal(:,1),FP.F_Signal(:,2), 'og', 'linewidth', 2)
plot(FP.D_Signal(:,1),FP.D_Signal(:,2),'oc','linewidth', 2)
str = ['Processed 470',filename];
    title(str) 
    xlabel('Time (minutes)')
    ylabel('z-score')
%%
% Find non-Interaction Epochs
  FP.B_Signal(:,1)=timeCorrected(:,1);
  FP.B_Signal(:,2)=FP.zdFF2(:,1);
  FP.B_SignalB(:,1)=timeCorrected(:,1);
  FP.B_SignalB(:,2)=FP.zdFF2(:,1);
  
  for i=1:length(FP.F_SegsMerge)
   start=find(FP.B_SignalB(:,1) > FP.F_SegsMerge(i,1)-0.0007 & FP.B_SignalB(:,1) < FP.F_SegsMerge(i,1)+0.0007);
   start=start(1,1);
   stop=find(FP.B_SignalB(:,1) > FP.F_SegsMerge(i,2)-0.0007 & FP.B_SignalB(:,1) < FP.F_SegsMerge(i,2)+0.0007);
   stop=stop(1,1);
   FP.B_Signal(start:stop,1)=zeros; 
   FP.B_Signal(start:stop,2)=zeros;  
  end
  
   for i=1:length(FP.D_SegsMerge)
   start=find(FP.B_SignalB(:,1) > FP.D_SegsMerge(i,1)-0.0007 & FP.B_SignalB(:,1) < FP.D_SegsMerge(i,1)+0.0007);
   start=start(1,1);
   stop=find(FP.B_SignalB(:,1) > FP.D_SegsMerge(i,2)-0.0007 & FP.B_SignalB(:,1) < FP.D_SegsMerge(i,2)+0.0007);
   stop=stop(1,1);
   FP.B_Signal(start:stop,1)=zeros; 
   FP.B_Signal(start:stop,2)=zeros;  
    end
     
  FP.B_Signal2(:,1)=nonzeros(FP.B_Signal(:,1));
  FP.B_Signal2(:,2)=nonzeros(FP.B_Signal(:,2));
    
  figure
  plot(FP.B_Signal2(:,1),FP.B_Signal2(:,2))
  xlim([0 5.2]);
  hold on
  plot(FP.F_Signal(:,1),FP.F_Signal(:,2), 'og', 'linewidth', 2)
  plot(FP.D_Signal(:,1),FP.D_Signal(:,2),'oc','linewidth', 2)
  
       
    
%%
%Perform analysis of all segments
Parameters.F_meanzdFF=mean(FP.F_Signal(:,2)); %Overall mean zdFF of F segments
Parameters.D_meanzdFF=mean(FP.D_Signal(:,2)); %Overall mean zdFF of D segments
Parameters.B_meanzdFF=mean(FP.B_Signal2(:,2)); %Overall mean zdFF of B (non) segments
Parameters.F_time=length(FP.F_Signal)*0.48; %F Time
Parameters.D_time=length(FP.D_Signal)*0.48; %D Time
Parameters.B_time=length(FP.B_Signal2)*0.48; %B Time

[D.y,D.x,D.w,D.p] = findpeaks(smooth(FP.D_Signal(:,2),10),'MinPeakProminence',0.04); % specify minimum prominence
FP.D_Peaks(:,1)=FP.D_Signal((D.x),1);
FP.D_Peaks(:,2)=D.y;
    
scatter(FP.D_Peaks(:,1),FP.D_Peaks(:,2),'k')

[F.y,F.x,F.w,F.p] = findpeaks(smooth(FP.F_Signal(:,2),10),'MinPeakProminence',0.04); % specify minimum prominence
FP.F_Peaks(:,1)=FP.F_Signal((F.x),1);
FP.F_Peaks(:,2)=F.y;
    
scatter(FP.F_Peaks(:,1),FP.F_Peaks(:,2),'k')

[B.y,B.x,B.w,B.p] = findpeaks(smooth(FP.B_Signal2(:,2),10),'MinPeakProminence',0.04); % specify minimum prominence
FP.B_Peaks(:,1)=FP.B_Signal2((B.x),1);
FP.B_Peaks(:,2)=B.y;
    
scatter(FP.B_Peaks(:,1),FP.B_Peaks(:,2),'b')

Parameters.F_meanpeakzdFF = mean(F.y);
Parameters.F_meanWidth = mean(F.w)*0.048;
Parameters.F_meanProm = mean(F.p);
Parameters.D_meanpeakzdFF = mean(D.y);
Parameters.D_meanWidth = mean(D.w)*0.048;
Parameters.D_meanProm = mean(D.p);
Parameters.B_meanpeakzdFF = mean(B.y);
Parameters.B_meanWidth = mean(B.w)*0.048;
Parameters.B_meanProm = mean(B.p);

Parameters.All(1,1) = Parameters.D_meanzdFF;
Parameters.All(2,1) = Parameters.D_meanpeakzdFF;
Parameters.All(3,1) = Parameters.D_meanWidth;
Parameters.All(4,1) = Parameters.D_meanProm;
Parameters.All(5,1) = Parameters.D_time;

Parameters.All(1,2) = Parameters.F_meanzdFF;
Parameters.All(2,2) = Parameters.F_meanpeakzdFF;
Parameters.All(3,2) = Parameters.F_meanWidth;
Parameters.All(4,2) = Parameters.F_meanProm;
Parameters.All(5,2) = Parameters.F_time;

Parameters.All(1,3) = Parameters.B_meanzdFF;
Parameters.All(2,3) = Parameters.B_meanpeakzdFF;
Parameters.All(3,3) = Parameters.B_meanWidth;
Parameters.All(4,3) = Parameters.B_meanProm;
Parameters.All(5,3) = Parameters.B_time;

annotation('textbox',[.906 .6 .1 .2],'String','D','EdgeColor','none','Color','blue')
annotation('textbox',[.906 .55 .1 .2],'String',Parameters.D_meanzdFF,'EdgeColor','none')
annotation('textbox',[.906 .5 .1 .2],'String',Parameters.D_meanpeakzdFF,'EdgeColor','none')
annotation('textbox',[.906 .45 .1 .2],'String',Parameters.D_meanWidth,'EdgeColor','none')
annotation('textbox',[.906 .4 .1 .2],'String',Parameters.D_meanProm,'EdgeColor','none')
annotation('textbox',[.906 .35 .1 .2],'String',Parameters.D_time,'EdgeColor','none')

annotation('textbox',[.94 .6 .1 .2],'String','F','EdgeColor','none','Color','green')
annotation('textbox',[.94 .55 .1 .2],'String',Parameters.F_meanzdFF,'EdgeColor','none')
annotation('textbox',[.94 .5 .1 .2],'String',Parameters.F_meanpeakzdFF,'EdgeColor','none')
annotation('textbox',[.94 .45 .1 .2],'String',Parameters.F_meanWidth,'EdgeColor','none')
annotation('textbox',[.94 .4 .1 .2],'String',Parameters.F_meanProm,'EdgeColor','none')
annotation('textbox',[.94 .35 .1 .2],'String',Parameters.F_time,'EdgeColor','none')

annotation('textbox',[.975 .6 .1 .2],'String','Non','EdgeColor','none','Color','black')
annotation('textbox',[.975 .55 .1 .2],'String',Parameters.B_meanzdFF,'EdgeColor','none')
annotation('textbox',[.975 .5 .1 .2],'String',Parameters.B_meanpeakzdFF,'EdgeColor','none')
annotation('textbox',[.975 .45 .1 .2],'String',Parameters.B_meanWidth,'EdgeColor','none')
annotation('textbox',[.975 .4 .1 .2],'String',Parameters.B_meanProm,'EdgeColor','none')
annotation('textbox',[.975 .35 .1 .2],'String',Parameters.B_time,'EdgeColor','none')

      set(gcf, 'Position', [100, 400, 1800, 250])
      xlim([0 5.2]);
      fig5=gcf;
     fig5str='/fig5.jpg';
     fig5bstr='/fig5.fig';
     filefig5str=strcat(filename,fig5str);
     filefig5bstr=strcat(filename,fig5bstr);
    saveas(fig5,filefig5str,'jpeg'); %save figure
    saveas(fig5,filefig5bstr,'fig');

    clear fig5
    clear fig5str
    clear fig5bstr
    clear filefig5str
    clear filefig5bstr
    
       
figure
subplot(2,1,1)
plot(FP.deinterleave(:,1),FP.deinterleave(:,2),'m')
xlim([0 5.2]);
 str = ['470-(',filename];
    title(str) 
    xlabel('Time (minutes)')
    ylabel('Mean Pixel Value')
  
    
subplot(2,1,2)    
plot(timeCorrected(:,1),FP.zdFF2(:,1),'k');
xlim([0 5.2]);
hold on
plot(FP.F_Signal(:,1),FP.F_Signal(:,2),'og','linewidth', 2)
plot(FP.D_Signal(:,1),FP.D_Signal(:,2),'oc','linewidth', 2)

scatter(FP.D_Peaks(:,1),FP.D_Peaks(:,2),'k')
scatter(FP.F_Peaks(:,1),FP.F_Peaks(:,2),'k')
scatter(FP.B_Peaks(:,1),FP.B_Peaks(:,2),'b')
    
str = ['Processed 470',filename];
    title(str) 
    xlabel('Time (minutes)')
    ylabel('z-score')
    
    set(gcf, 'Position', [400, 300, 1000, 600])
     fig6=gcf;
     fig6str='/fig6.jpg';
     fig6bstr='/fig6.fig';
     filefig6str=strcat(filename,fig6str);
     filefig6bstr=strcat(filename,fig6bstr);
    saveas(fig6,filefig6str,'jpeg'); %save figure
    saveas(fig6,filefig6bstr,'fig');

    clear fig6
    clear fig6str
    clear fig6bstr
    clear filefig6str
    clear filefig6bstr
    
%%
%do analysis for each segment

%find mean zDFFs (D Segs)
    for i=1:length(FP.D_SegsMerge)
    number=num2str(i);
      name=strcat('FP.dzdFF2_',number);
      name2=strcat('Parameters.DSegs(',number);
      name2=strcat(name2,',1)');
      str=[name2,'=mean(',name,'(:,1));'];
      eval(str);
      clear str;
    end
    
%find mean zDFFs (F Segs)
        for i=1:length(FP.F_SegsMerge)
    number=num2str(i);
      name=strcat('FP.fzdFF2_',number);
      name2=strcat('Parameters.FSegs(',number);
      name2=strcat(name2,',1)');
      str=[name2,'=mean(',name,'(:,1));'];
      eval(str);
      clear str;
        end
        
                %do peaks analysis (D Segs)
    for i=1:length(FP.D_SegsMerge)
    number=num2str(i);
      name=strcat('FP.dzdFF2_',number);
      name2=strcat('Parameters.DSegs(',number);
      str=['[Dt.y,Dt.x,Dt.w,Dt.p] = findpeaks(smooth(',name,'(:,1),10),"MinPeakProminence",0.05);'];
      eval(str);
      clear str;
      Parameters.DSegs(i,2)=mean(Dt.w)*0.048;
      Parameters.DSegs(i,3)=mean(Dt.p);
      clear Dt
    end
        
        %do peaks analysis (F Segs)
    for i=1:length(FP.F_SegsMerge)
    number=num2str(i);
      name=strcat('FP.fzdFF2_',number);
      name2=strcat('Parameters.FSegs(',number);
      str=['[Ft.y,Ft.x,Ft.w,Ft.p] = findpeaks(smooth(',name,'(:,1),10),"MinPeakProminence",0.05);'];
      eval(str);
      clear str;
      Parameters.FSegs(i,2)=mean(Ft.w)*0.048;
      Parameters.FSegs(i,3)=mean(Ft.p);
      clear Ft
    end
    
    figure
    set(gcf, 'Position', [400, 300, 1000, 600])
    subplot(2,3,1)
    plot(Parameters.DSegs(:,1),'c');
   
    title('D Mean zdFF') 
    xlabel('Segment')
    ylabel('zdFF')
    
    subplot(2,3,2)
    plot(Parameters.DSegs(:,2),'c');
    
    title('D Mean Peak Width') 
    xlabel('Segment')
    ylabel('time (s)')
    
    subplot(2,3,3)
    plot(Parameters.DSegs(:,3),'c');

    title('D Mean Peak Prom') 
    xlabel('Segment')
    ylabel('Prominence (z)')
    
        subplot(2,3,4)
    plot(Parameters.FSegs(:,1),'g');

    title('F Mean zdFF') 
    xlabel('Segment')
    ylabel('zdFF')
    
    subplot(2,3,5)
    plot(Parameters.FSegs(:,2),'g');

    title('F Mean Peak Width') 
    xlabel('Segment')
    ylabel('time (s)')
    
    subplot(2,3,6)
    plot(Parameters.FSegs(:,3),'g');

    title('F Mean Peak Prom') 
    xlabel('Segment')
    ylabel('Prominence (z)')
    

     fig7=gcf;
     fig7str='/fig7.jpg';
     fig7bstr='/fig7.fig';
     filefig7str=strcat(filename,fig7str);
     filefig7bstr=strcat(filename,fig7bstr);
    saveas(fig7,filefig7str,'jpeg'); %save figure
    saveas(fig7,filefig7bstr,'fig');

    clear fig7
    clear fig7str
    clear fig7bstr
    clear filefig7str
    clear filefig7bstr
    
     %save(filename)
      
      sized=size(Parameters.DSegs);
      sizef=size(Parameters.FSegs);
      
      
      if sized > sizef
          sizew=sized;
      else
          sizew=sizef;
      end
      
      Parameters.All((end+1):(sizew(1,1)),1:3)=zeros;
      Parameters.All(1:(sizew(1,1)),4:9)=zeros;
      Parameters.All(1:(sized(1,1)),4:6)=Parameters.DSegs(:,:);
      Parameters.All(1:(sizef(1,1)),7:9)=Parameters.FSegs(:,:);
      
      % write results to csv
% 
  paramstr='/parameters.csv';
  paramstr2=strcat(filename,paramstr);
  csvwrite(paramstr2, Parameters.All);
  
  sizecheckd=size(FP.D_Signal(:,1));
  sizecheckf=size(FP.F_Signal(:,1));
  sizecheckb=size(FP.B_Signal2(:,1));
  sizecheckO=size(FP.zdFF2(:,1));
  sizecheckA=sizecheckd+sizecheckf+sizecheckb;
  sizecheckA(1,1) == sizecheckO(1,1);
  
  
      
      %%
    %time-locked analysis, various baseline correction approaches were
    %tested
    
     for i=1:length(FP.D_SegsMerge)
     
      FP.dSignal=find(timeCorrected(:,1) > FP.D_SegsMerge(i,1)-0.0007 & timeCorrected (:,1) < FP.D_SegsMerge(i,1)+0.0007);
      FP.D_Segments(:,i)=FP.zdFF2(FP.dSignal(1,1)-120:FP.dSignal(1,1)+120,1);
      clear FP.dSignal
     end
        
      figure
      for i=1:size(FP.D_Segments,2)
          hold on
          plot(FP.D_Segments(:,i));
      end
            
      for i=1:size(FP.D_Segments,2)
          FP.D_SegmentsDFF(:,i) = (FP.D_Segments(:,i)/mean(FP.D_Segments(1:120,i)));
      end
      
            figure
      for i=1:size(FP.D_SegmentsDFF,2)
          hold on
          plot(FP.D_SegmentsDFF(:,i));
      end
      
      for i=1:size(FP.D_Segments,2)
          FP.D_SegmentsBL(:,i) = (FP.D_Segments(:,i) + (0-(FP.D_Segments(120,i))));
      end
      
     
      FP.D_SegsBLmean(:,1)=mean(FP.D_SegmentsBL,2);
      FP.D_SegsMean(:,1)=mean(FP.D_Segments,2);

      
      for i=1:size(FP.D_Segments,2)
          FP.D_SegmentsBL2(:,i) = (FP.D_Segments(:,i) + (0-min(FP.D_Segments(:,i))));
      end
      
            FP.D_SegsBL2mean(:,1)=mean(FP.D_SegmentsBL2,2);
            
            
                          for i=1:size(FP.D_Segments,2)
          FP.D_SegmentsBL3(:,i) = (FP.D_Segments(:,i) + (0-(FP.D_Segments(50,i))));
                          end
                
                          
            
                         for i=1:size(FP.D_Segments,2) %this is the method for baselining used
          FP.D_SegmentsBL3B(:,i) = (FP.D_Segments(:,i) + (0-(mean(FP.D_Segments(30:80,i)))));
                          end
                
                   
            FP.D_SegsBL3mean(:,1)=mean(FP.D_SegmentsBL3,2);
            
            FP.D_SegsBL3Bmean(:,1)=mean(FP.D_SegmentsBL3B,2);
            
 figure
      for i=1:size(FP.D_SegmentsDFF,2)
          hold on
          plot(FP.D_SegmentsBL(:,i));
      end
                 
 figure
      for i=1:size(FP.D_SegmentsDFF,2)
          hold on
          plot(FP.D_SegmentsBL2(:,i));
      end
      %%
      
            for i=1:length(FP.F_SegsMerge)
     
      FP.fSignal=find(timeCorrected(:,1) > FP.F_SegsMerge(i,1)-0.0007 & timeCorrected (:,1) < FP.F_SegsMerge(i,1)+0.0007);
      FP.F_Segments(:,i)=FP.zdFF2(FP.fSignal(1,1)-120:FP.fSignal(1,1)+120,1);
      clear FP.fSignal
            end
            
                     figure
      for i=1:size(FP.F_Segments,2)
          hold on
          plot(FP.F_Segments(:,i));
      end
      
       for i=1:size(FP.F_Segments,2)
          FP.F_SegmentsDFF(:,i) = (FP.F_Segments(:,i)/mean(FP.F_Segments(1:120,i)));
      end
      
            figure
      for i=1:size(FP.F_SegmentsDFF,2)
          hold on
          plot(FP.F_SegmentsDFF(:,i));
      end
      
                for i=1:size(FP.F_Segments,2)
          FP.F_SegmentsBL(:,i) = (FP.F_Segments(:,i) + (0-(FP.F_Segments(120,i))));
                end
          
          for i=1:size(FP.F_Segments,2)
          FP.F_SegmentsBL2(:,i) = (FP.F_Segments(:,i) + (0-min(FP.F_Segments(:,i))));
          end
      
                          for i=1:size(FP.F_Segments,2)
          FP.F_SegmentsBL3(:,i) = (FP.F_Segments(:,i) + (0-(FP.F_Segments(50,i))));
                          end
                   
                
            
                         for i=1:size(FP.F_Segments,2)
          FP.F_SegmentsBL3B(:,i) = (FP.F_Segments(:,i) + (0-(mean(FP.F_Segments(30:80,i)))));
                          end
                
                   
          
                      figure
      for i=1:size(FP.F_SegmentsDFF,2)
          hold on
          plot(FP.F_SegmentsBL(:,i));
      end
            
      
      FP.F_SegsBLmean(:,1)=mean(FP.F_SegmentsBL,2);
      FP.F_SegsBL2mean(:,1)=mean(FP.F_SegmentsBL2,2);
      FP.F_SegsMean(:,1)=mean(FP.F_Segments,2);
      
            FP.F_SegsBL3mean(:,1)=mean(FP.F_SegmentsBL3,2);
            
            FP.F_SegsBL3Bmean(:,1)=mean(FP.F_SegmentsBL3B,2);
                 
                      figure
      for i=1:size(FP.F_SegmentsDFF,2)
          hold on
          plot(FP.F_SegmentsBL2(:,i));
      end
      
      figure
      hold on
      plot(FP.D_SegsBLmean)
      plot(FP.F_SegsBLmean)
      
      
      figure
      hold on
      plot(FP.D_SegsBL2mean)
      plot(FP.F_SegsBL2mean)
      
      
      figure
      hold on
      plot(FP.D_SegsBL3mean)
      plot(FP.F_SegsBL3mean)
      
      %%
      FP.zdFF3=FP.zdFF2 * std(signal) + median(raw_signal);
      
    
     for i=1:length(FP.D_SegsMerge)
     
      FP.dSignal=find(timeCorrected(:,1) > FP.D_SegsMerge(i,1)-0.0007 & timeCorrected (:,1) < FP.D_SegsMerge(i,1)+0.0007);
      FP.D_SegmentsALT(:,i)=FP.zdFF3(FP.dSignal(1,1)-120:FP.dSignal(1,1)+120,1);
      clear FP.dSignal
     end
      
      
            for i=1:length(FP.F_SegsMerge)
     
      FP.fSignal=find(timeCorrected(:,1) > FP.F_SegsMerge(i,1)-0.0007 & timeCorrected (:,1) < FP.F_SegsMerge(i,1)+0.0007);
      FP.F_SegmentsALT(:,i)=FP.zdFF3(FP.fSignal(1,1)-120:FP.fSignal(1,1)+120,1);
      clear FP.fSignal
            end
            
            
      for i=1:size(FP.D_SegmentsALT,2)
          FP.D_SegmentsDFFALT(:,i) = ((FP.D_SegmentsALT(:,i) - mean(FP.D_SegmentsALT(30:80,i)))/mean(FP.D_SegmentsALT(30:80,i)));
      end
      
      
      for i=1:size(FP.F_SegmentsALT,2)
          FP.F_SegmentsDFFALT(:,i) = ((FP.F_SegmentsALT(:,i) - mean(FP.F_SegmentsALT(30:80,i)))/mean(FP.F_SegmentsALT(30:80,i)));
      end
      
      
            FP.D_SegsBL4mean(:,1)=mean(FP.D_SegmentsDFFALT,2);
            FP.F_SegsBL4mean(:,1)=mean(FP.F_SegmentsDFFALT,2);
            
            figure
            hold on
            plot(FP.D_SegsBL4mean)
            plot(FP.F_SegsBL4mean)
                 
    %%
            
    figure
    subplot(1,5,1)
      hold on
      plot(FP.D_SegsBLmean)
      plot(FP.F_SegsBLmean)
    subplot(1,5,2)
      hold on
      plot(FP.D_SegsBL2mean)
      plot(FP.F_SegsBL2mean)
          subplot(1,5,3)
      hold on
      plot(FP.D_SegsBL3mean)
      plot(FP.F_SegsBL3mean)
            subplot(1,5,4)
      hold on
      plot(FP.D_SegsBL3Bmean)
      plot(FP.F_SegsBL3Bmean)
                subplot(1,5,5)
      hold on
      plot(FP.D_SegsBL4mean)
      plot(FP.F_SegsBL4mean)
      
                 set(gcf, 'Position', [200, 300, 1700, 320])
     figTemp=gcf;
     figTempStr='/TimeLocked.jpg';
     figTempBStr='/TimeLocked.fig';
     filefigTempStr=strcat(filename,figTempStr);
     filefigTempBStr=strcat(filename,figTempBStr);
    saveas(figTemp,filefigTempStr,'jpeg'); %save figure
    saveas(figTemp,filefigTempBStr,'fig');
    
    clear figTemp
    
    
    %%
    
    FP.TimeLocked(:,1)=FP.D_SegsBLmean;
    FP.TimeLocked(:,2)=FP.F_SegsBLmean;
    
    FP.TimeLocked(:,3)=FP.D_SegsBL2mean;
    FP.TimeLocked(:,4)=FP.F_SegsBL2mean;
    
    FP.TimeLocked(:,5)=FP.D_SegsBL3mean;
    FP.TimeLocked(:,6)=FP.F_SegsBL3mean;
     
    FP.TimeLocked(:,7)=FP.D_SegsBL3Bmean; %this is the baseline method used
    FP.TimeLocked(:,8)=FP.F_SegsBL3Bmean;
    
    FP.TimeLocked(:,9)=FP.D_SegsBL4mean;
    FP.TimeLocked(:,10)=FP.F_SegsBL4mean;
    
      save(filename)
      % 
  timestr='/timelocked.csv';
  timestr2=strcat(filename,timestr);
  csvwrite(timestr2, FP.TimeLocked);
  

  