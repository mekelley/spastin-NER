%outputs of interest: 1) set 1, set 2 are cells containing tracks associated to each pole with
%col1: dist, col2: time, col3: associated pole ID. Each cell contains a
%diff track. 2) DistPole2Pole_Time contains the pole to pole distance
%versus time
clear all; close all;
pxl=1;
px=1;


[DataFile,DataPath]=uigetfile('*.csv','Select tracks file', 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_cenCENP\raw_crops\SPZ\Tracking');
    Data=xlsread(fullfile(DataPath,DataFile)); %input
%   Data=xlsread('C:\Users\linac\Documents\MATLAB\MEKSpastin\kinetocTracks\202016_RPECen2x_0p1DMSO_10min_2s_backsub_TM_SPOTS.csv');

newName=split(DataFile,'.nd2');
%newName=split(DataFile,'.csv');

seriesLong=split(DataFile,'series ');
series=split(seriesLong(2),')');

%coreName=newName(1);
coreName=newName(1);

Tracks_cell=[Data(:,4), Data(:,5), Data(:,6), Data(:,8), Data(:,2)]; %x,y,z, frame,id, new version trackmate
deltaT=30; %input
min_length=8;
matchTime=8;
max_length=10000;

prefactor=0.10;

res_filt=[];
one=[];two=[];three=[];four=[]; five=[];

for i=0:max(Tracks_cell(:,5))+1
     
    
    [target, ~]=find(Tracks_cell(:,5)==i);
    
    
   if length(target)>=min_length && length(target)<=max_length  
        
    one=cat(1,one,Tracks_cell(target,1));
    two=cat(1,two,Tracks_cell(target,2));
    three=cat(1,three,Tracks_cell(target,3));
    four=cat(1,four,Tracks_cell(target,4));
    five=cat(1,five,Tracks_cell(target,5));
    
    %% one, two three , four, five are filtered
   else
   clear target
    end
end
pxl=1;
% FiltTracks are the filtered tracks
FiltTracks(:,5)=five; %track ID
FiltTracks(:,4)=deltaT*four; %time (s)
FiltTracks(:,1)=pxl*one;%x (um)
FiltTracks(:,2)=pxl*two;%y (um)
FiltTracks(:,3)=pxl*three;
  clear one two three four
  
  MaxTime=max(FiltTracks(:,4));
  ViewTime=prefactor*MaxTime; %change the prefactor if you can't find the poles

count=1;
for i=0:max(FiltTracks(:,5))+1
     
    [target, ~] = find(FiltTracks(:,5)==i);
    
    
    tracks{count,1}(:,1)=FiltTracks(target,1); % x
    tracks{count,1}(:,2)=FiltTracks(target,2); % y
    tracks{count,1}(:,3)=FiltTracks(target,3); % z
        tracks{count,1}(:,4)=FiltTracks(target,4); % time

    
    count=count+1;
    end

tracks_v=tracks(~cellfun('isempty',tracks));
Tracks2=tracks_v';
 count=0;
for i=1:length(Tracks2)
   
   
    if Tracks2{i}(end,4)<=MaxTime % silly condition, but keep for now. might be useful later
         [target, ~]=find(Tracks2{i}(:,4)<=ViewTime);
        count=1+count;
        TracksToView{count}=Tracks2{i}(target, :);
    else
        
    end
end

 %get velocity 
%  Tracks2=Tracks(idx);
 for t=1:length(Tracks2)
     for r=1:(size(Tracks2{1,t},1)-1)
         %instantaneous velocity
        v{t}(r,1)=px.*(sqrt((Tracks2{1,t}(r+1,1)-Tracks2{1,t}(r,1)).^2+(Tracks2{1,t}(r+1,2)-Tracks2{1,t}(r,2)).^2+(Tracks2{1,t}(r+1,3)-Tracks2{1,t}(r,3)).^2))./(Tracks2{1,t}(r+1,4)-Tracks2{1,t}(r,4));
     end
   
     v_mean(t)=mean(v{t}(:,:));
     
 end
 
 close all;
 for i=1:length(Tracks2)

% scatter3(Tracks2{1,i}(:,1),Tracks2{1,i}(:,2), Tracks2{1,i}(:,3), 'ro')
hold on
plot3(TracksToView{1,i}(:,1),TracksToView{1,i}(:,2), TracksToView{1,i}(:,3))
 end
 title('Choose poles')
%  view(3)
 %choose the centrosome
[CentX, CentY]=getpts;

close all
for i=1:length(Tracks2)

% scatter3(Tracks2{1,i}(:,1),Tracks2{1,i}(:,2), Tracks2{1,i}(:,3), 'ro')
hold on
plot3(TracksToView{1,i}(:,1),TracksToView{1,i}(:,2), TracksToView{1,i}(:,3))
end
 title('Select tracks to remove. If none, hit enter')
%  view(3)
 %choose the centrosome
[removeX, removeY]=getpts;



if ~isempty(removeX)
     for j=1:length(removeX)
for i=1:length(TracksToView)
    
    [IDremove(i), D_remove(i)]=knnsearchEXT([removeX(j), removeY(j)],TracksToView{1,i}(:,1:2));
  
end
  id=find(D_remove==min(D_remove));
  id;
  TracksToView{1,id}=zeros(0,4);
  Tracks2{1,id}=zeros(0,4);
  clear id
     end
end
close

 for i=1:length(Tracks2)

% scatter3(Tracks2{1,i}(:,1),Tracks2{1,i}(:,2), Tracks2{1,i}(:,3), 'ro')
hold on
plot3(TracksToView{1,i}(:,1),TracksToView{1,i}(:,2), TracksToView{1,i}(:,3))
 end
 title('removal check')
 
 
%identify tracks closest to user selection 1 (poles)
for ii=1:length(Tracks2)
    
    [IDX_1(ii), D_1(ii)]=knnsearchEXT([CentX(1), CentY(1)],Tracks2{1,ii}(:,1:2));
    [IDX_2(ii), D_2(ii)]=knnsearchEXT([CentX(2), CentY(2)],Tracks2{1,ii}(:,1:2)); 
   
    %verify this is also in the short list of tracks
    [IDX_1v(ii), D_1v(ii)]=knnsearchEXT([CentX(1), CentY(1)],TracksToView{1,ii}(:,1:2));
    [IDX_2v(ii), D_2v(ii)]=knnsearchEXT([CentX(2), CentY(2)],TracksToView{1,ii}(:,1:2)); 
    
end

overlap1_binary=ismember(IDX_1, IDX_1v);
overlap2_binary=ismember(IDX_2, IDX_2v);

d_1=overlap1_binary.*D_1;
d_2=overlap2_binary.*D_2;

d_1(find(d_1==0))=NaN;
d_2(find(d_2==0))=NaN;

 CentTrack_id1 = find(d_1==min(d_1)); % first centrosome track
 CentTrack_id2 = find(d_2==min(d_2)); % second centrosome track


 %find centrosome that each track is closest to
 for t=1:length(Tracks2)
     if ~isempty(Tracks2{1,t})
%      timer = Tracks2{1,t}(end,4);
timer = Tracks2{1,t}(matchTime,4);
     temp1 = find(Tracks2{1,CentTrack_id1}(:,4)==timer);
     temp2 = find(Tracks2{1,CentTrack_id2}(:,4)==timer);
    
     if isempty(temp1) || isempty(temp2) 
         temp1=size(Tracks2{1,CentTrack_id1}, 2);
         temp2=size(Tracks2{1,CentTrack_id1}, 2);
     end
      
%      if isempty(temp2)==1
%           temp2=size(Tracks2{1,CentTrack_id1}, 2);
%      end
%      
%      IDXchoice(t)=knnsearchEXT(Tracks2{1,t}(end,1:2),[Tracks2{1,CentTrack_id1}(temp1,1:2);Tracks2{1,CentTrack_id2}(temp2,1:2)]);
          IDXchoice(t)=knnsearchEXT(Tracks2{1,t}(matchTime,1:2),[Tracks2{1,CentTrack_id1}(temp1,1:2);Tracks2{1,CentTrack_id2}(temp2,1:2)]);

     if IDXchoice(t)==1
         realchoice(t)= CentTrack_id1;
     else if IDXchoice(t)==2
             realchoice(t)= CentTrack_id2;
         end
     end
     else 
         if isempty(Tracks2{1,t})
             realchoice(t)=NaN;
         end
     end
     clear temp 1 temp2 timer
 end
 
 poleIDs=unique(realchoice);


%consider only times when both poles are measured
 [idxPole1, ~]=find(ismember(Tracks2{poleIDs(1)}(:,4),Tracks2{poleIDs(2)}(:,4))); %time indices for pole 1
 [idxPole2, ~]=find(ismember(Tracks2{poleIDs(2)}(:,4),Tracks2{poleIDs(1)}(:,4))); %time indices for pole 2, not the lengths of these idxPoles are identical
  
deltaX = Tracks2{poleIDs(1)}(idxPole1,1)-Tracks2{poleIDs(2)}(idxPole2,1);
deltaY =  Tracks2{poleIDs(1)}(idxPole1,2)-Tracks2{poleIDs(2)}(idxPole2,2); 
deltaZ =  Tracks2{poleIDs(1)}(idxPole1,3)-Tracks2{poleIDs(2)}(idxPole2,3);

dist=sqrt(deltaX.^2+deltaY.^2+deltaZ.^2);
 DistPole2Pole_Time= [dist, Tracks2{poleIDs(1)}(idxPole1,4)];
  close all;
  figure;hold on;
  plot(DistPole2Pole_Time(:,2), DistPole2Pole_Time(:,1), 'k-'); %Output of interest where Col 2 is time and Col1 is pole-to-pole distance
  xlabel('Time [s]')
  ylabel('Pole-to-pole distance [\mum]')


% for t=1:length(Tracks2{poleIDs(OtherPole)})
%     DistPole2Pole_T(t)

%for each kinetochore track:
% 1) align by time
% 2) calculate kinetochore to pole distance
 for m=1:length(Tracks2)
     m;
%      kTrack=size(Tracks2{m},1);
%      CentTrack=size(Tracks2{realchoice(m)},1);
%      if kTrack<CentTrack 
%          idxLong=find(ismember(Tracks2{realchoice(m)}(:,4),Tracks2{m}(:,4)));
if ~isnan(realchoice(m))
         [idxLong1, ~] = find(ismember(Tracks2{realchoice(m)}(:,4),Tracks2{m}(:,4))); %this ensures dist measurement is taken only when pole track exists
         [idxLong2, ~] = find(ismember(Tracks2{m}(:,4),Tracks2{realchoice(m)}(:,4))); %same as comment in line above 
         
        


         for r=1:length(idxLong1)%(size(Tracks2{1,m},1))
             if m~=realchoice(m)
             DeltaK2Pole{m}=(Tracks2{m}(idxLong2,1:3)-Tracks2{realchoice(m)}(idxLong1,1:3)); %first input is kineto, second is pole 
             DistK2Pole_Time{m}(r,:)= [sqrt(DeltaK2Pole{m}(r,1).^2+DeltaK2Pole{m}(r,2).^2+DeltaK2Pole{m}(r,3).^2), Tracks2{m}(idxLong2(r),4), realchoice(m)];
             end
         end
         
         
else 
end
        
         
%      else if kTrack>=CentTrack
%          idxLong=find(ismember(Tracks2{m}(:,4),Tracks2{realchoice(m)}(:,4)));
%      for r=1:length(idxLong);%size(Tracks2{realchoice(m)},1)
%              DeltaK2Pole{m}=(Tracks2{realchoice(m)}(:,1:3)-Tracks2{m}(idxLong,1:3));
%              DistK2Pole_Time{m}(r,:)= [sqrt(DeltaK2Pole{m}(r,1).^2+DeltaK2Pole{m}(r,2).^2+DeltaK2Pole{m}(r,3).^2), idxLong(r)];
%          end
%      end
%      clear kTrack; clear CentTrack; clear idxLong
 
 end
      Dist2pole_T=DistK2Pole_Time(~cellfun('isempty',DistK2Pole_Time)); %Dist2pole_T is output of interest where the col 2 is time and col1 is dist bw k and pole

 count1=0;
 count2=0;
 for i =1: length(Dist2pole_T)
     
     if Dist2pole_T{i}(1,3) == poleIDs(1)
     count1=count1+1;
     set1{count1}=Dist2pole_T{i};
     else 
         count2=count2+1;
          set2{count2}=Dist2pole_T{i};
     end
 end
     figure;
 for i=1:length(set1)
hold on
plot(set1{i}(:,2),set1{i}(:,1), 'r-')
end
hold on
for i=1:length(set2)
hold on
plot(set2{i}(:,2),set2{i}(:,1), 'b-')
end
 
     figure;
     
     for i=1:length(Dist2pole_T)
         
         hold on
         plot(Dist2pole_T{i}(:,2), Dist2pole_T{i}(:,1), '-') 
     end

xlabel('Time [s]')
ylabel('Kinetochore to pole distance [\mum]')
    
 
 MeanSpeed=mean(v_mean); %instantaneous mean speed
 SpeedEr=std(v_mean); %std deviation of distribution of inst. mean speed
 
cd(DataPath);
%saveFILE=strcat(char(coreName),'.mat');
saveFILE=strcat(char(coreName),'_',char(series(1)),'.mat');
save(saveFILE);


% 
% 
% for i=1:length(Dist2pole_T)
%     if isempty(Dist2pole_T{1,i})~=1
%         if Dist2pole_T{1,i}(1,1)>=10
%          i
%         end
%     end
% end
% 
%       Dist2pole_T=Dist2pole_T(~cellfun('isempty',Dist2pole_T)); %Dist2pole_T is output of interest where the col 2 is time and col1 is dist bw k and pole
% 
% for i=1:length(DistK2Pole_Time)
%     if isempty(DistK2Pole_Time{1,i})~=1
%         if DistK2Pole_Time{1,i}(1,1)>=10
%          i
%         end
%     end
% end