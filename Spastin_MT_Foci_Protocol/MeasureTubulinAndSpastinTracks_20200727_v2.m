clear all;
close all;
[TracksFile,TracksPath]=uigetfile('*.csv','Select Tracks File', 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Tracks');
csvFile=fullfile(TracksPath,TracksFile);
Data=xlsread(csvFile);
Tracks_cell=[Data(:,4), Data(:,5), Data(:,8), Data(:,2)]; %x,y,frame,id, new version trackmate

cd('C:\Users\kelle\Documents\MATLAB\saveastiff_4.4')
[C1file,C1filePath]=uigetfile('*.tif','Select TUB C1 File', 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Sub70\Split');
TubFile=fullfile(C1filePath,C1file);
Tub=loadtiff(TubFile);
[C2file,C2filePath]=uigetfile('*.tif','Select SPAST C2 File', 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Sub70\Split');
SpasFile=fullfile(C2filePath,C2file);
Spas=loadtiff(SpasFile);

delay=15; %frame delay between 'real' tubulin frames
deltaT=1; %time delay in s
pxl= (1./9.0909);%0.1148;
r1=4; %radius over which to integrate intensity of spastin and tubulin signal
r2=8; %r2-r1 is the ring-shaped area over which spastin background is calculated
r3=r1;


TubFrames=1:delay:size(Tub,3);
if TubFrames(end)>size(Tub,3)
    TubFrames(end)=[];
else
end

%filter tracks
min_length=5;
max_length=1000;
res_filt=[];
one=[];two=[];three=[];four=[];

for i=0:max(Tracks_cell(:,4))+1
     
    
    target=find(Tracks_cell(:,4)==i);
    
    
   if length(target)>=min_length && length(target)<=max_length  
        
    one=cat(1,one,Tracks_cell(target,1));
    two=cat(1,two,Tracks_cell(target,2));
    three=cat(1,three,Tracks_cell(target,3));
    four=cat(1,four,Tracks_cell(target,4));
    
    %% one, two three and four are filtered
   else
   clear target
    end
end

% FiltTracks are the filtered tracks
FiltTracks(:,1)=four; %track ID
FiltTracks(:,2)=deltaT*three; %time (s)
FiltTracks(:,3)=1+one./pxl;%x (um)
FiltTracks(:,4)=1+two./pxl;%y (um)

  clear one two three four

count=1;
for i=0:max(FiltTracks(:,1))+1
     
    target=find(FiltTracks(:,1)==i);
    
    if length(target)>min_length
    
    tracks{count,1}(:,1)=1+FiltTracks(target,2); % frame
    tracks{count,1}(:,2)=FiltTracks(target,3); % x in pixel
     tracks{count,1}(:,3)=FiltTracks(target,4); % y in pixels
    tracks{count, 1}(:,4)=FiltTracks(target,1); %trackmate ID
    count=count+1;
    else
    end
end

%Tubulin MIP for background definition
    TubMIP1=uint16(zeros(size(Tub(:,:,1))));
    for j=1:round(size(Tub,3)/2)
    TubMIP1=max(TubMIP1, uint16(Tub(:,:,j)));
    end
    
      TubMIP2=uint16(zeros(size(Tub(:,:,1))));
    for j=round(size(Tub,3)/2):size(Tub,3)
    TubMIP2=max(TubMIP2, uint16(Tub(:,:,j)));
    end
    
    imshow(TubMIP1,[0, .7*max(max(TubMIP1))], 'initialmagnification', 300)
hold on;

%user defined background area
BackgroundPoly=impoly;
BackgroundPos=BackgroundPoly.getPosition;
 BackgroundPos1(:,1)=[BackgroundPos(:,1); BackgroundPos(1,1)];
 BackgroundPos1(:,2)=[BackgroundPos(:,2); BackgroundPos(1,2)];
close;
 mask1(:,:)=poly2mask(BackgroundPos1(:,1),BackgroundPos1(:,2), size(TubMIP1, 1), size(TubMIP1,2));
 
 clear BackgroundPoly BackgroundPos BackgroundPos1 
 
 imshow(TubMIP2,[0, .7*max(max(TubMIP2))], 'initialmagnification', 300)
hold on;
BackgroundPoly=impoly;
BackgroundPos=BackgroundPoly.getPosition;
 BackgroundPos1(:,1)=[BackgroundPos(:,1); BackgroundPos(1,1)];
 BackgroundPos1(:,2)=[BackgroundPos(:,2); BackgroundPos(1,2)];
close;
 mask2(:,:)=poly2mask(BackgroundPos1(:,1),BackgroundPos1(:,2), size(TubMIP2, 1), size(TubMIP2,2));


 for i=1:round(size(Tub,3)/2)
     BkTub(:,i)=max(max(immultiply(Tub(:,:,i), mask1))); %background for each frame
 end
 
 for i=round(size(Tub,3)/2):size(Tub,3)
     BkTub(:,i)=max(max(immultiply(Tub(:,:,i), mask2))); 
 end
 close all;
 plot(BkTub);
 
% define circles around spastin
q=0:0.01:2*pi; % angles to generate circles
for i=1:length(tracks)
    for j=1:size(tracks{i},1)
    SpasCircles{i}{j}=[(r1*cos(q)+tracks{i}(j,2))' (r1*sin(q)+tracks{i}(j,3))']; 
    BackgroundCircles{i}{j}=[(r2*cos(q)+tracks{i}(j,2))' (r2*sin(q)+tracks{i}(j,3))']; 
    end
end

% define circles around spastin for tubulin
q=0:0.01:2*pi; % angles to generate circles
for i=1:length(tracks)
    for j=1:size(tracks{i},1)
    SpasCircles_forTubulin{i}{j}=[(r3*cos(q)+tracks{i}(j,2))' (r3*sin(q)+tracks{i}(j,3))']; 
    end
end

%make binary masks (with circles above) and calculate the integrated
%intensity inside and local background 

for i=1:length(tracks)
    for j=1:size(tracks{i},1)
%         i
%         j
        tempInnerMask = poly2mask(SpasCircles{i}{j}(:,1),SpasCircles{i}{j}(:,2), size(Tub(:,:,1), 1), size(Tub(:,:,1), 2) );
        tempOuterMask = poly2mask(BackgroundCircles{i}{j}(:,1),BackgroundCircles{i}{j}(:,2), size(Tub(:,:,1), 1), size(Tub(:,:,1), 2) );
        tempTaurMask = imsubtract(tempOuterMask, tempInnerMask);
        
         tempInnerMaskTub = poly2mask(SpasCircles_forTubulin{i}{j}(:,1),SpasCircles_forTubulin{i}{j}(:,2), size(Tub(:,:,1), 1), size(Tub(:,:,1), 2) );

        clear tempOuterMask
        
        productImSpas = immultiply(double(Spas(:,:,tracks{i}(j,1))), double(tempInnerMask)); %product of real Spastin image and mask
        productImSpas_bk = immultiply(double(Spas(:,:,tracks{i}(j,1))), double(tempTaurMask)); %local background spastin
        productImTub =  immultiply(double(Tub(:,:,tracks{i}(j,1))), double(tempInnerMaskTub)); %product of real tubulin image and mask
%         productImTub_bk = immultiply(double(Tub(:,:,tracks{i}(j,1))), double(tempTaurMask)); %local background tubulin
        
        % for spastin:integrated int, max int, min int, max_bk, min_bk, frame
        SpasSignals{i}(j,1:6) = [sum(sum(productImSpas)), max(max(productImSpas)), min(min(productImSpas)), max(max(productImSpas_bk)), min(min(productImSpas_bk)), tracks{i}(j,1)];
               
        % for tubulin: integrated int, max int, min int, max int- bk, min int - bk, frame
        TubSignals{i}(j,1:6) = [sum(sum(productImTub)), max(max(productImTub)),min(min(productImTub)),(max(max(productImTub))-BkTub(tracks{i}(j,1))),(min(min(productImTub))-BkTub(tracks{i}(j,1))), tracks{i}(j,1)];
        
        clear tempInnerMask tempTaurMask  productImSpas  productImSpas_bk productImTub  productImTub_bk tempInnerMaskTub
    end
end

for i=1:length(TubSignals)
    
   if sum(TubSignals{i}(:,4))<=0
       BinTub=0;
   elseif sum(TubSignals{i}(:,4))>0
       BinTub=1;
   end
   TubPres(i)=BinTub; %presence of tubulin, 1=present at some point, 0 = never
   AllTubSignal_bkCorr(i)=sum(TubSignals{i}(:,4));
   clear BinTub
end
cutoff=prctile(AllTubSignal_bkCorr, 50); 
cutoff

PercColoc=length(find(TubPres==0))./length(find(TubPres==1)); %this is the percentage of punctae that colocalize with Tubulin (at some point in spot's lifetime)

idHighTub=find(AllTubSignal_bkCorr>=cutoff);
idLowTub=find(AllTubSignal_bkCorr<cutoff);

idTub=find(TubPres==1);
idNoTub=find(TubPres==0);

for i=1:length(idTub)
    meanTotal_wTub(i)=mean(SpasSignals{idTub(i)}(:,1));
    meanMax_wTub(i)=mean(SpasSignals{idTub(i)}(:,2));
    max_wTub(i)=max(SpasSignals{idTub(i)}(:,2));
    Tau_wTub(i)=size(SpasSignals{idTub(i)},1);
end

for i=1:length(idHighTub) % 'high' intensity of tubulin
    meanTotal_HTub(i)=mean(SpasSignals{idHighTub(i)}(:,1));
    meanMax_HTub(i)=mean(SpasSignals{idHighTub(i)}(:,2));
    max_HTub(i)=max(SpasSignals{idHighTub(i)}(:,2));
    Tau_HTub(i)=size(SpasSignals{idHighTub(i)},1);
end

for i=1:length(idNoTub)
    meanTotal_woTub(i) = mean(SpasSignals{idNoTub(i)}(:,1));
    meanMax_woTub(i) = mean(SpasSignals{idNoTub(i)}(:,2));
    max_woTub(i) = max(SpasSignals{idNoTub(i)}(:,2));
    Tau_woTub(i) = size(SpasSignals{idNoTub(i)},1);
end

for i=1:length(idLowTub) % 'Low' intensity of tubulin
    meanTotal_LTub(i)=mean(SpasSignals{idLowTub(i)}(:,1));
    meanMax_LTub(i)=mean(SpasSignals{idLowTub(i)}(:,2));
    max_LTub(i)=max(SpasSignals{idLowTub(i)}(:,2));
    Tau_LTub(i)=size(SpasSignals{idLowTub(i)},1);
end

%the following three cells will be populated with matrices containing: 
%col 1: real frame, col 2: spastin signal, col 3: tubulin signal with gap
%filling, col 4: 'real' tubulin signal, col 5: track id from trackmate file
NonCoreSignals={};
MidzoneCoreSignals={};
KinCoreSignals={};

counter1=0;
counter2=0;
counter3=0;


%% plotting and user-selection  (to do)   
close all
for i=1:length(SpasSignals)
 
    figure;
    yyaxis left
    plot(SpasSignals{i}(:,end),SpasSignals{i}(:,2), 'bo-' );
    hold on  
   
    yyaxis right
    plot(TubSignals{i}(:,end),TubSignals{i}(:,4), 'ro-' );
    hold on
    if isempty(intersect(TubSignals{i}(:,end), TubFrames))==0
        [value,pos]=intersect(TubSignals{i}(:,end), TubFrames);
        hold on
        yyaxis right
        plot(TubSignals{i}(pos,end),TubSignals{i}(pos,4), 'g*' )
    else
    end
          
legend('Spastin signal', 'Background corrected Tubulin')
            xlabel('Time [frame ]')
ylabel('Max Signal')

clear TempTimeProj

%time projections and tracks

  TempTimeProj=uint16(zeros(size(Tub(:,:,1))));
    for j=1:size(tracks{i},1)
    TempTimeProj=max( TempTimeProj, Tub(:,:,SpasSignals{i}(j,end)));
    end
    
figure;    
    imshow(TempTimeProj, [], 'initialmagnification', 300); hold on; 
    for t=1:size(tracks{i},1)
        plot(SpasCircles{i}{t}(:,1),SpasCircles{i}{t}(:,2), 'm-' )
    end
%     hold on
%     plot(tracks{i}(:,2),tracks{i}(:,3), 'm-' )



% %plot localizations over individual frames
% for t=1:size(tracks{i},1)
% figure;
% imshowpair(Spas(:,:,SpasSignals{i}(t,end)),Tub(:,:,SpasSignals{i}(t,end)), 'montage' );
% hold on
% plot(SpasCircles{i}{t}(:,1),SpasCircles{i}{t}(:,2), 'w.' )
% hold on
% plot(size(Tub(:,:,1), 1)+SpasCircles{i}{t}(:,1),SpasCircles{i}{t}(:,2),'w.' ); %displaying only first frame at which track appears
pause

choice = questdlg('Is this envelope associated?', ...
    'Choice', ...
    'No','Yes!','No');


% Handle response
switch choice
    case 'No'
    case 'Yes!'
        choice2 = questdlg('Designate region', ...
              'Choice', ...
    'Non-core','Kinetochore', 'Midzone', 'Non-core');

    switch choice2

    case 'Non-core'
        counter1=counter1+1;
        
        RealTubITemp=zeros(length(SpasSignals{i}(:,end)), 1); %initialize a matrix to input real intensity values corresponding to the real tubulin signal
        clear pos value
         if isempty(intersect(TubSignals{i}(:,end), TubFrames))==0
        [value,pos]=intersect(TubSignals{i}(:,end), TubFrames);
        for ii=1:length(pos)
            RealTubITemp(pos(ii))=TubSignals{i}(pos(ii),4);
        end
    else
    end
        
        NonCoreSignals{counter1} = [SpasSignals{i}(:,end), SpasSignals{i}(:,2),TubSignals{i}(:,4), RealTubITemp, tracks{i}(:,4)];
        
        clear RealTubITemp
        
    case 'Kinetochore'
        counter2=counter2+1;
        
        
        RealTubITemp=zeros(length(SpasSignals{i}(:,end)), 1); %initialize a matrix to input real intensity values corresponding to the real tubulin signal
        clear pos value
         if isempty(intersect(TubSignals{i}(:,end), TubFrames))==0
        [value,pos]=intersect(TubSignals{i}(:,end), TubFrames);
        for ii=1:length(pos)
            RealTubITemp(pos(ii))=TubSignals{i}(pos(ii),4);
        end
    else
    end
        
        KinCoreSignals{counter2} = [SpasSignals{i}(:,end), SpasSignals{i}(:,2),TubSignals{i}(:,4),   RealTubITemp, tracks{i}(:,4)];
                clear RealTubITemp

        
    case 'Midzone'
        counter3=counter3+1;
        
            RealTubITemp=zeros(length(SpasSignals{i}(:,end)), 1); %initialize a matrix to input real intensity values corresponding to the real tubulin signal
        clear pos value
         if isempty(intersect(TubSignals{i}(:,end), TubFrames))==0
        [value,pos]=intersect(TubSignals{i}(:,end), TubFrames);
        for ii=1:length(pos)
            RealTubITemp(pos(ii))=TubSignals{i}(pos(ii),4);
        end
    else
        end
    
            MidzoneCoreSignals{counter3} = [SpasSignals{i}(:,end), SpasSignals{i}(:,2),TubSignals{i}(:,4),   RealTubITemp, tracks{i}(:,4)];
                    clear RealTubITemp


    end
end
pause
close all
end

%if you quit selection, you can run from here (FYI)
     close all;
     MidzoneCore_Tub=[];
         MidzoneCore_noTub=[];
     KinCore_Tub=[];
          KinCore_noTub=[];
     Noncore=[];
 clear i
 for i=1:length(MidzoneCoreSignals)
     if sum(MidzoneCoreSignals{i}(:,4))>0
        MidzoneCore_Tub=[MidzoneCore_Tub; mean(MidzoneCoreSignals{i}(:,2))  max(MidzoneCoreSignals{i}(:,4)) size(MidzoneCoreSignals{i}(:,2), 1) MidzoneCoreSignals{i}(1,5)] ; %mean spas signal, max tub, lifetime, trackmate ID
     elseif sum(MidzoneCoreSignals{i}(:,4))==0
        MidzoneCore_noTub=[MidzoneCore_noTub; mean(MidzoneCoreSignals{i}(:,2))  max(MidzoneCoreSignals{i}(:,4)) size(MidzoneCoreSignals{i}(:,2), 1) MidzoneCoreSignals{i}(1,5)] ; %mean spas signal, max tub, lifetime, trackmate ID
     end
 end 
 
 for i=1:length(KinCoreSignals)
       if sum(KinCoreSignals{i}(:,4))>0
     KinCore_Tub=[KinCore_Tub; mean(KinCoreSignals{i}(:,2))  max(KinCoreSignals{i}(:,4)) size(KinCoreSignals{i}(:,2), 1) KinCoreSignals{i}(1,5)] ;   %mean spas signal, max tub, lifetime, trackmate ID
       elseif sum(KinCoreSignals{i}(:,4))==0
                 KinCore_noTub=[KinCore_noTub; mean(KinCoreSignals{i}(:,2))  max(KinCoreSignals{i}(:,4)) size(KinCoreSignals{i}(:,2), 1) KinCoreSignals{i}(1,5)] ; %mean spas signal, max tub, lifetime, trackmate ID
         end
 end
 
 
 for i=1:length(NonCoreSignals)
     Noncore=[Noncore; mean(NonCoreSignals{i}(:,2))  max(NonCoreSignals{i}(:,4)) size(NonCoreSignals{i}(:,2), 1) NonCoreSignals{i}(1,5)];    %mean spas signal, max tub, lifetime, trackmate ID
 end
 
%plot lifetimes
histogram(Noncore(:,end))
hold on; histogram(MidzoneCore(:,end-1)); hold on; histogram(KinCore(:,end-1))
title('Lifetimes')

cd(Folder);
saveFILE=strcat(char(saveName), '.mat');
save(saveFILE);
     
     