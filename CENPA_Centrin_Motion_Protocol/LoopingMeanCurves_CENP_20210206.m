clearvars -except stotalExps DtotalExps
%clear all
%close all

% Specify the folder where the files live.
theFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_cenCENP\raw_crops\SPZ\Tracking';
%StheFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_cenCENP\raw_crops\SPZ\Tracking';

testTitle = '20210223_testCENP';
cd('D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_cenCENP\raw_crops');

d1spz=['1 ',char(181),'M spastazoline'];
d2p5spz=['2.5 ',char(181),'M spastazoline'];
d5spz=['5 ',char(181),'M spastazoline'];
d10spz=['10 ',char(181),'M spastazoline'];
wtspz=strcat('WT',{' '},d1spz);
N386Cspz=strcat('N386C',{' '},d1spz);
micron=[char(181),'m'];

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(theFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', theFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name
% pattern.c
fPattern = fullfile(theFolder, '*.mat'); % Change to whatever pattern you need.
matFiles = dir(fPattern);

for n = 1:length(matFiles) 
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    baseName = matFiles(n).name;
    fullName = fullfile(theFolder, baseName);
    [filepath,name,ext] = fileparts(fullName);
    matName=name;
    %spName = split(name,'_');
    
    load(fullfile(filepath, strcat(name,ext)));
    close all;
    
    k=1;
    for i=1:size(DistK2Pole_Time,2)
        if isempty(DistK2Pole_Time{1,i}==1)
        else
            AllData{k,1}=DistK2Pole_Time{1,i};
            k=k+1;
        end
    end
    
%  for i=1:size(DistK2Pole_Time,2)
%         k=DistK2Pole_Time(1,i);
%         AllData2{k,1}(n,1:3)=DistK2Pole_Time(:,i);
%     end
    
        BINsize=30; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,2));
          ShortestTime(t)=min(AllData{t}(:,2));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,2)>=TimeStamp(i) & AllData{ii}(:,2)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,1))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=mean(BinIntInt{m});
              stdIntInt(m)=std(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,1}=matName;
    totalExps{n,2}=TimeStamp;
    totalExps{n,3}=meanIntInt;
    totalExps{n,4}=stdIntInt;
    totalExps{n,5}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp
end





%%%SUM OF SUMMARY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for DMSO to be in black and spastazoline in magenta
figure;
for n=1:length(matFiles)
    if mod(n,2)~=0 
        
        
        hold on
%         title('Centroid Distance')
        shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'bl', 1);
        %shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'm', 1);
        xlabel('Time (sec)') 
        ylabel('Distance (microns)')
        %axis([0 100 0 70])
    end   
end

figure;
for n=1:length(matFiles)
    if mod(n,2)~=0 
        
        
        hold on
%         title('Centroid Distance')
        %shadedErrorBar(DtotalExps{n,2},DtotalExps{n,3},DtotalExps{n,4}, 'bl', 1);
        shadedErrorBar(stotalExps{n,2},stotalExps{n,3},stotalExps{n,4}, 'm', 1);
        xlabel('Time (sec)') 
        ylabel('Distance (microns)')
        %axis([0 100 0 70])
    end   
end

figure;
for n=1:length(DtotalExps)
    if mod(n,2)~=0 
        
        
        hold on
%         title('Centroid Distance')
        shadedErrorBar(DtotalExps{n,2},DtotalExps{n,3},DtotalExps{n,4}, 'bl', 1);

        xlabel('Time (sec)') 
        ylabel('Distance (microns)')
        %axis([0 100 0 70])
    end   
end
for n=1:length(DtotalExps)
    if mod(n,2)~=0 
        
        
        hold on
%         title('Centroid Distance')
        shadedErrorBar(stotalExps{n,2},stotalExps{n,3},stotalExps{n,4}, 'm', 1);
        %axis([0 100 0 70])
    end   
end



